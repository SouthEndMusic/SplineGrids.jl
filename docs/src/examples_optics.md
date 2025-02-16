# Lens geometry optimization example

In this section we optimize the surface of a lens to yield a target irradiance distribution on a screen. This process is called _caustic design_. Below a schematic of the simulated setup is shown.

```@example
using CairoMakie # hide
using Makie.GeometryBasics # hide

fig = Figure() # hide
ax = Axis(fig[1, 1]; xlabel = "z axis") # hide

# Incoming beam # hide
n_rays = 10 # hide
arrows!(ax, # hide
    fill(-2.0, n_rays), # hide
    range(-1, 1, length = n_rays), # hide
    ones(n_rays),  # hide
    zeros(n_rays), # hide
    label = "Incoming parallel beam" # hide
) # hide

# Lens first surface # hide
n_points = 100  # hide
y_lens = range(-1, 1, length = n_points) # hide
x_lens_1 = fill(-0.5, n_points) # hide
lines!(ax, # hide
    x_lens_1,  # hide
    y_lens, # hide
    label = "Fixed lens surface" # hide
) # hide

# Lens second surface # hide
x_lens_2 = sin.(10y_lens) / 25  # hide
lines!(ax, # hide
    x_lens_2, # hide
    y_lens, # hide
    label = "Lens surface to optimize" # hide
) # hide

# Lens fill # hide
poly!( # hide
    ax, # hide
    [ # hide
        Point2f.(x_lens_1, y_lens)..., # hide
        Point2f.(x_lens_2[end:-1:1], y_lens[end:-1:1])... # hide
    ], # hide
    color = (:blue, 0.2) # hide
) # hide

# Detector screen # hide
lines!( # hide
    ax, # hide
    [2.0, 2.0], # hide
    [-3.0, 3.0], # hide
    label = "Detector screen" # hide
) # hide

axislegend(ax, position = :lt) # hide
hidexdecorations!(ax, label = false) # hide
hideydecorations!(ax) # hide

fig # hide
```

## Differentiable ray tracing kernel

Here we define a kernel which computes for each sample point on the surface:

  - at which direction a ray leaves that point by Snell's law;
  - where that ray intersects the detector screen;
  - What the contributions of that ray to the pixel values are.

The contribution of the ray is smeared out over multiple pixels in a smooth way to make the rendering differentiable.

```@example tutorial
using KernelAbstractions
using Atomix

# See the included file for the ray tracing implementation
include("ray_tracing_kernel.jl")
```

## Calling the ray tracing kernel

Here we define a function which computes the input for the ray tracing kernel from a spline grid and then calls the kernel.

```@example tutorial
function trace_rays!(render, control_points_flat, p)::Nothing
    (; spline_grid, u, ∂₁u, ∂₂u) = p

    copyto!(spline_grid.control_points, control_points_flat)
    evaluate!(spline_grid.control_points) # Only relevant for LocallyRefinedControlPoints

    evaluate!(spline_grid; eval = u)
    evaluate!(spline_grid; eval = ∂₁u, derivative_order = (1, 0))
    evaluate!(spline_grid; eval = ∂₂u, derivative_order = (0, 1))

    render .= 0.0
    backend = get_backend(u)

    ray_tracing_kernel(backend)(
        render,
        u,
        ∂₁u,
        ∂₂u,
        spline_grid.spline_dimensions[1].sample_points,
        spline_grid.spline_dimensions[2].sample_points,
        p.scene_params,
        ndrange = size(u)
    )
    synchronize(backend)
    return nothing
end
```

## Tracing the first rays

Let's define a flat spline surface and trace some rays. We expect to see a projection of the square lens onto the screen, as all rays travel parallel to the z-axis.

```@example tutorial
using SplineGrids
using Plots
using CairoMakie # hide

n_control_points = (20, 20)
degree = (2, 2)
n_sample_points = (1000, 1000) # Determines grid of sampled rays
dim_out = 1
extent = (-1.0, 1.0) # Lens extent in both x and y direction

spline_dimensions = SplineDimension.(
    n_control_points, degree, n_sample_points; max_derivative_order = 1, extent)
spline_grid = SplineGrid(spline_dimensions, dim_out)
spline_grid.control_points .= 0
fig_total = Figure(size = (1600, 900)) # hide
plot_basis!(fig_total, spline_grid; j = 2, title = "Unrefined spline basis") # hide

# Static scene parameters
scene_params = (;
    r = 1.4f0,
    z_screen = 4.0f0,
    screen_size = (7.5f0, 7.5f0),
    ray_kernel_size = (5, 5),
    screen_res = (250, 250)
)

p_render = (;
    spline_grid,
    u = similar(spline_grid.eval),
    ∂₁u = similar(spline_grid.eval),
    ∂₂u = similar(spline_grid.eval),
    scene_params
)

render = zeros(Float32, scene_params.screen_res)

trace_rays!(render, vec(spline_grid.control_points), p_render)

Plots.heatmap(render, aspect_ratio = :equal, cmap = :grays,
    title = "Initial distribution with flat lens")
```

## Defining the target distribution

We define a normalized target distribution, which we will compare to normalized renders.

```@example tutorial
target = [exp(-(x .^ 2 + y .^ 2)^2)
          for
          x in range(-2, 2,
    length = scene_params.screen_res[1]
),
y in range(-2, 2,
    length = scene_params.screen_res[2]
)]

# Normalize target
target ./= sum(target)

ax_target = Axis(fig_total[2, 1]; title = "Target illumination", aspect = 1) # hide
CairoMakie.heatmap!(ax_target, target, colormap = :grays) # hide
Plots.heatmap(target, aspect_ratio = :equal, cmap = :grays, title = "Target distribution")
```

## The loss function

```@example tutorial
using Distances

function image_loss(control_points_flat, target, render, p_render)
    trace_rays!(render, control_points_flat, p_render)
    # normalize render
    render ./= sum(render)
    Euclidean()(render, target)
end

image_loss(
    vec(spline_grid.control_points),
    target,
    render,
    p_render
)
```

## Gradients w.r.t. control points

We can now compute the gradient of the loss function with respect to the control points. Let's have a look at it.

```@example tutorial
using Enzyme

G = make_zero(vec(spline_grid.control_points))
drender = make_zero(render)
dp_render = make_zero(p_render)

autodiff(
    Enzyme.Reverse,
    image_loss,
    Active,
    Duplicated(vec(spline_grid.control_points), G),
    Const(target),
    DuplicatedNoNeed(render, drender),
    DuplicatedNoNeed(p_render, dp_render)
)

M = maximum(abs.(G))
Plots.heatmap(
    reshape(G, n_control_points), aspect_ratio = :equal, cmap = :bluesreds, clims = (-M, M),
    title = "The gradient of the loss function\nw.r.t. the control point grid")
```

## Optimizing the surface

```@example tutorial
using Optimization
using OptimizationOptimJL: Adam

function image_loss_grad!(G, control_points_flat, meta_p)::Nothing
    make_zero!(G)
    (; render_duplicated, p_render_duplicated) = meta_p
    (; spline_grid, u, ∂₁u, ∂₂u) = p_render_duplicated.dval
    make_zero!(render_duplicated.dval)
    make_zero!(spline_grid)
    make_zero!(u)
    make_zero!(∂₁u)
    make_zero!(∂₂u)
    autodiff(
        Enzyme.Reverse,
        image_loss,
        Active,
        Duplicated(control_points_flat, G),
        Const(meta_p.target),
        meta_p.render_duplicated,
        meta_p.p_render_duplicated
    )
    return nothing
end

meta_p = (;
    target,
    render_duplicated = DuplicatedNoNeed(render, drender),
    p_render_duplicated = DuplicatedNoNeed(p_render, dp_render)
)

optimization_function = OptimizationFunction(
    (control_points_flat, p) -> image_loss(
        control_points_flat,
        target,
        render,
        p_render
    ),
    grad = image_loss_grad!
)

prob = OptimizationProblem(
    optimization_function,
    vec(spline_grid.control_points),
    meta_p
)

loss_values = Float32[]
function loss_callback(state, args...)
    (; objective) = state
    push!(loss_values, objective)
    println("Loss (iteration #$(length(loss_values))) = $objective")
    # Stop the optimization if the loss doesn't reduce
    false
end

sol = solve(prob, Adam(; alpha = 1f-3, epsilon=1f-8, beta_mean=0.9f0, beta_var=0.999f0); maxiters = 20, callback = loss_callback)
```

## Viewing the optimization result

The loss reduces over the iterations as follows:

```@example tutorial
Plots.plot(loss_values, yscale = :log10, xlabel = "Iterations", title = "Loss")
```

The final render looks like this:

```@example tutorial
trace_rays!(render, sol.u, p_render)
ax_render = Axis(
    fig_total[2, 2]; aspect = 1, title = "Illumination after first optimization") # hide
CairoMakie.heatmap!(ax_render, render; colormap = :grays) # hide
#err =  # hide
Plots.heatmap(render, aspect_ratio = :equal, cmap = :grays, title = "Final render")
```

And the lens surface looks like this:

```@example tutorial
spline_grid.control_points .= reshape(sol.u, size(spline_grid.control_points))
evaluate!(spline_grid)
M = maximum(abs.(spline_grid.eval))
ax_lens = Axis(fig_total[4, 2]; aspect = 1, title = "Lens surface after first optimization") # hide
CairoMakie.heatmap!(ax_lens, spline_grid.eval[:, :, 1]) # hide
err = render ./ sum(render) - target # hide
ax_err = Axis(fig_total[3, 2]; aspect = 1, title = "Error after first optimization")
M = maximum(abs.(err))
CairoMakie.heatmap!(ax_err, err, colorrange = (-M, M), colormap = :redblue)
Plots.plot(spline_grid; plot_knots = false, aspect_ratio = :equal,
    cmap = :viridis, title = "Final lens surface")
```

## Locally refining the basis of the lens surface

One of the neat things we can do with this setup is look at all sorts of gradients. We are most interested in the gradient of the loss with respect to the partial derivatives of the surface, since those are the most important for the rendering result. In particular, we look at the sum of the absolute values of these gradients. This shows which regions of the lens surface the loss is most sensitive to, and thus where the surface might need more degrees of freedom.

```@example tutorial
function loss_from_grid(render, u, ∂₁u, ∂₂u, spline_grid, target, p_render)
    backend = get_backend(u)
    ray_tracing_kernel(backend)(
        render,
        u,
        ∂₁u,
        ∂₂u,
        spline_grid.spline_dimensions[1].sample_points,
        spline_grid.spline_dimensions[2].sample_points,
        p_render.scene_params,
        ndrange = size(u)
    )
    synchronize(backend)
    Euclidean()(render, target)
end

for val in values(meta_p.p_render_duplicated.dval)
    val isa Union{Array, SplineGrid} && make_zero!(val)
end
autodiff(
    Enzyme.Reverse,
    loss_from_grid,
    Active,
    meta_p.render_duplicated,
    Duplicated(p_render.u, meta_p.p_render_duplicated.dval.u),
    Duplicated(p_render.∂₁u, meta_p.p_render_duplicated.dval.∂₁u),
    Duplicated(p_render.∂₂u, meta_p.p_render_duplicated.dval.∂₂u),
    Duplicated(spline_grid, meta_p.p_render_duplicated.dval.spline_grid),
    Const(target),
    Const(p_render)
)

sensitivity = abs.(meta_p.p_render_duplicated.dval.∂₁u) +
              abs.(meta_p.p_render_duplicated.dval.∂₂u)

Plots.heatmap(
    sensitivity[:, :, 1], aspect_ratio = :equal, title = "Loss sensitivity of lens surface")
```

With the local refinement functionality of `SplineGrids.jl`, we can use this sensitivity data to refine the spline basis of the lens surface in those regions where the sensitivity to the loss is highest:

```@example tutorial
using CairoMakie

fig = Figure()
plot_basis!(fig, spline_grid, title = "Unrefined basis")

spline_grid = add_default_local_refinement(spline_grid)
error_informed_local_refinement!(spline_grid, sensitivity)
deactivate_overwritten_control_points!(spline_grid.control_points)
plot_basis!(fig_total, spline_grid; j = 3, title = "Refined spline basis") # hide
plot_basis!(fig, spline_grid; j = 2, title = "Refined basis")
fig
```

Now let's see whether the optimization can do better.

```@example tutorial
using ConstructionBase
control_points = zeros(Float32, get_n_control_points(spline_grid), 1)
copyto!(control_points, spline_grid.control_points)

p_render = setproperties(p_render; spline_grid)
meta_p = setproperties(
    meta_p; p_render_duplicated = DuplicatedNoNeed(p_render, make_zero(p_render)))

prob = OptimizationProblem(
    optimization_function,
    control_points,
    meta_p
)

sol = solve(prob, Adam(; alpha = 1f-3, epsilon=1f-8, beta_mean=0.9f0, beta_var=0.999f0); maxiters = 100, callback = loss_callback)
```

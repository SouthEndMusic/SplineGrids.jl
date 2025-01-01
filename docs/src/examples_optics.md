# Lens geometry optimization example

In this section we optimize the surface of a lens to yield a target irradiance distribution on a screen. This process is called _caustic design_. Below a schematic of the simulated setup is shown.

```@example
using CairoMakie # hide
using Makie.GeometryBasics # hide

fig = Figure() # hide
ax = Axis(fig[1,1]; xlabel = "z axis") # hide

# Incoming beam # hide
n_rays = 10 # hide
arrows!(ax, # hide
    fill(-2.0, n_rays), # hide
    range(-1,1, length = n_rays), # hide
    ones(n_rays),  # hide
    zeros(n_rays), # hide
    label = "Incoming parallel beam" # hide
) # hide

# Lens first surface # hide
n_points = 100  # hide
y_lens = range(-1,1, length = n_points) # hide
x_lens_1 = fill(-0.5, n_points) # hide
lines!(ax, # hide
    x_lens_1,  # hide
    y_lens, # hide
    label = "Fixed lens surface" # hide
) # hide

# Lens second surface # hide
x_lens_2 = sin.(10y_lens)/25  # hide
lines!(ax, # hide
    x_lens_2, # hide
    y_lens, # hide
    label = "Lens surface to optimize", # hide
) # hide

# Lens fill # hide
poly!( # hide
    ax, # hide
    [ # hide
        Point2f.(x_lens_1, y_lens)..., # hide
        Point2f.(x_lens_2[end:-1:1], y_lens[end:-1:1])... # hide
    ], # hide
    color = (:blue, 0.2), # hide
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

# Kernel function for computing the contribution of a ray intersection to the
# pixels close to the intersection.
function F(x, x0, w)
    if x < x0 - w
        -one(x) / 2
    elseif x > x0 + w
        one(x) / 2
    else
        x_transformed = (x - x0) / w
        (sin(π * x_transformed) / π + x_transformed) / 2
    end
end

@kernel function ray_tracing_kernel(
    render,
    @Const(u),
    @Const(∂₁u),
    @Const(∂₂u),
    @Const(x),
    @Const(y),
    r,
    z_screen,
    screen_size,
    ray_kernel_size,
)

    I = @index(Global, Cartesian)

    u_I = u[I]
    ∂₁u_I = ∂₁u[I]
    ∂₂u_I = ∂₂u[I]

    x_I = x[I[1]]
    y_I = y[I[2]]

    # x direction tangent vector: (1, 0, ∂₁u[I])
    # y direction tangent vector: (0, 1, ∂₂u[I])
    # -cross product (surface normal): n = (∂₁u[I], ∂₂u[I], -1) / √(1 + ∂₁u[I]^2 + ∂₁u[I]^2)
    # light vector: ℓ = (0, 0, 1)
    # Snell's law (vector form):
    # c = -⟨n, ℓ⟩ = 1 / √(1 + ∂₁u[I]^2 + ∂₁u[I]^2)
    # v = r * (0, 0, 1) + (r * c - √(1 - r^2 * (1 - c^2))) * (∂₁u[I], ∂₂u[I], -1) / √(1 + ∂₁u[I]^2 + ∂₁u[I]^2)

    cross_product_norm = √(1 + ∂₁u_I^2 + ∂₂u_I^2)
    c = 1 / cross_product_norm
    sqrt_arg = 1 - r^2 * (1 - c^2)
    if sqrt_arg >= 0
        normal_vector_coef = (r * c - √(1 - r^2 * (1 - c^2))) / cross_product_norm

        # Refracted ray direction
        v_x = normal_vector_coef * ∂₁u_I
        v_y = normal_vector_coef * ∂₂u_I
        v_z = -normal_vector_coef + r

        # Refracted ray starting point: (x_I, y_I, u_I)
        t_screen_int = (z_screen - u_I) / v_z

        if t_screen_int >= 0

            # Screen intersection coordinates
            x_screen = x_I + t_screen_int * v_x
            y_screen = y_I + t_screen_int * v_y

            # Pixel size
            w_screen, h_screen = screen_size
            n_x, n_y = size(render)
            w_pixel = w_screen / n_x
            h_pixel = h_screen / n_y

            # Pixel intersection indices
            n_x, n_y = size(render)

            i = 1 + Int(floor((w_screen / 2 + x_screen) / w_pixel))
            j = 1 + Int(floor((h_screen / 2 + y_screen) / h_pixel))

            # Render contribution from this ray
            i_min = max(i - ray_kernel_size[1] - 1, 1)
            i_max = min(i + ray_kernel_size[1] + 1, n_x)
            j_min = max(j - ray_kernel_size[2] - 1, 1)
            j_max = min(j + ray_kernel_size[2] + 1, n_y)

            w_kernel = (ray_kernel_size[1] + 0.5) * w_pixel
            h_kernel = (ray_kernel_size[2] + 0.5) * h_pixel

            for i_ in i_min:i_max
                contribution_x =
                    F(-0.5w_screen + i_ * w_pixel, x_screen, w_kernel) -
                    F(-0.5w_screen + (i_ - 1) * w_pixel, x_screen, w_kernel)

                for j_ in j_min:j_max
                    contribution_y =
                        F(-0.5h_screen + j_ * h_pixel, y_screen, h_kernel) -
                        F(-0.5h_screen + (j_ - 1) * h_pixel, y_screen, h_kernel)

                    Atomix.@atomic render[i_, j_] += contribution_x * contribution_y
                end
            end
        end
    end
end
```

## Calling the ray tracing kernel

Here we define a function which computes the input for the ray tracing kernel from a spline grid and then calls the kernel.

```@example tutorial
function trace_rays!(render, control_points_flat, p)::Nothing
    (; spline_grid, u, ∂₁u, ∂₂u) = p

    control_points = reshape(control_points_flat, size(spline_grid.control_points))

    evaluate!(spline_grid; control_points, eval=u)
    evaluate!(spline_grid; control_points, eval=∂₁u, derivative_order=(1, 0))
    evaluate!(spline_grid; control_points, eval=∂₂u, derivative_order=(0, 1))

    render .= 0.0
    backend = get_backend(u)

    ray_tracing_kernel(backend)(
        render,
        u,
        ∂₁u,
        ∂₂u,
        spline_grid.spline_dimensions[1].sample_points,
        spline_grid.spline_dimensions[2].sample_points,
        p.r,
        p.z_screen,
        p.screen_size,
        p.ray_kernel_size,
        ndrange=size(u)
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

n_control_points = (50, 50)
degree = (2, 2)
n_sample_points = (300, 300) # Determines grid of sampled rays
dim_out = 1
extent = (-1.0, 1.0) # Lens extent in both x and y direction

spline_dimensions = SplineDimension.(n_control_points, degree, n_sample_points; max_derivative_order=1, extent)
spline_grid = SplineGrid(spline_dimensions, dim_out)
spline_grid.control_points .= 0

p_render = (;
    spline_grid,
    u=similar(spline_grid.eval),
    ∂₁u=similar(spline_grid.eval),
    ∂₂u=similar(spline_grid.eval),
    r=1.4,
    z_screen=5.0,
    screen_size=(4.0, 4.0),
    ray_kernel_size=(3, 3),
    screen_res=(250, 250)
)

render = zeros(Float32, p_render.screen_res)

trace_rays!(render, vec(spline_grid.control_points), p_render)

heatmap(render, aspect_ratio=:equal)
```

## Defining the target distribution

We define a normalized target distribution, which we will compare to normalized renders.

```@example tutorial
using LinearAlgebra

target = [
    exp(-(x.^2 + y.^2)^2) for
    x = range(-p_render.screen_size[1] / 2, p_render.screen_size[1] / 2, length=p_render.screen_res[1]),
    y = range(-p_render.screen_size[2] / 2, p_render.screen_size[2] / 2, length=p_render.screen_res[2])
]

normalize!(target)

heatmap(target, aspect_ratio=:equal)
```

## The loss function

```@example tutorial
using Distances

function image_loss(control_points_flat, target, render, p_render)

    trace_rays!(render, control_points_flat, p_render)
    normalize!(render)
    Euclidean()(render, target)
end

render = zeros(Float32, p_render.screen_res...)

image_loss(
    vec(spline_grid.control_points),
    target,
    render,
    p_render,
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
    Reverse,
    image_loss,
    Active,
    Duplicated(vec(spline_grid.control_points), G),
    Const(target),
    DuplicatedNoNeed(render, drender),
    DuplicatedNoNeed(p_render, dp_render),
)

heatmap(reshape(G, n_control_points), aspect_ratio=:equal)
```

## Optimizing the surface

```@example tutorial
using Optimization
using OptimizationOptimJL: BFGS

function image_loss_grad!(G, control_points_flat, meta_p)::Nothing
    make_zero!(G)
    make_zero!(meta_p.render_duplicated.dval)
    for val in values(meta_p.p_render_duplicated.dval)
        val isa Union{Array, SplineGrid} && make_zero!(val)
    end
    autodiff(
        Reverse,
        image_loss,
        Active,
        Duplicated(control_points_flat, G),
        Const(meta_p.target),
        meta_p.render_duplicated,
        meta_p.p_render_duplicated,
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
    meta_p,
)

sol = solve(prob, BFGS(); maxiters = 50)
```

## Viewing the optimization result

The final render looks like this:

```@example tutorial
trace_rays!(render, sol.u, p_render)
heatmap(render, aspect_ratio=:equal)
```

And the lens surface looks like this:

```@example tutorial
spline_grid.control_points .= reshape(sol.u, size(spline_grid.control_points))
evaluate!(spline_grid)
plot(spline_grid; plot_knots = false, aspect_ratio=:equal)
```

## A peek into upcoming features

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
        p_render.r,
        p_render.z_screen,
        p_render.screen_size,
        p_render.ray_kernel_size,
        ndrange=size(u)
    )
    synchronize(backend)
    Euclidean()(render, target)
end

for val in values(meta_p.p_render_duplicated.dval)
    val isa Union{Array, SplineGrid} && make_zero!(val)
end
autodiff(
    Reverse,
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

heatmap(
    abs.(meta_p.p_render_duplicated.dval.∂₁u[:,:,1]) +
    abs.(meta_p.p_render_duplicated.dval.∂₂u[:,:,1]), 
    aspect_ratio=:equal)
```
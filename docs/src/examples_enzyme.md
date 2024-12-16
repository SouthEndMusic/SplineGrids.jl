# Enzyme example

To demonstrate how Enzyme can be used to incorporate a spline grid into an optimization pipeline, we will fit a spline surface to the image below.

```@example tutorial
using Images
image = load(normpath(@__DIR__, "julia_logo.png"))
```

## Defining the spline grid

We define a spline grid with 2 input dimensions, 1 output dimension and a sample grid which matches the image resolution.

```@example tutorial
using SplineGrids

n_control_points = (40, 40)
degree = (2, 2)
image_array = Float64.(Gray.(image[end:-1:1, :]))'
n_sample_points = size(image_array)
dim_out = 1

spline_dimensions = SplineDimension.(n_control_points, degree, n_sample_points)
spline_grid = SplineGrid(spline_dimensions, dim_out)
spline_grid.control_points .= 0
spline_grid
```

## Gradient demonstration

Here we show what the gradient of the output surface with respect to one control point looks like. 

```@example tutorial
using Enzyme
using Plots

dspline_grid = deepcopy(spline_grid)
dspline_grid.control_points[5, 5] = 1

out = autodiff(Forward, Duplicated(spline_grid, dspline_grid)) do spline_grid
    evaluate!(spline_grid)
    spline_grid.eval
end

heatmap(out[1][:, :, 1])
```

## Fitting

We fit the spline surface to the image by minimizing a loss function.

```@example tutorial
using Optimization
using OptimizationOptimJL: BFGS

# Loss function
function image_loss(control_points_flat, spline_grid)
    control_points = reshape(control_points_flat, size(spline_grid.control_points))
    evaluate!(spline_grid; control_points)
    out = zero(eltype(image_array))
    for I in eachindex(image_array)
        out += (spline_grid.eval[I] - image_array[I])^2
    end
    return sqrt(out / length(image_array))
end

# Loss function gradient
function image_loss_grad!(G, control_points_flat, p)::Nothing
    make_zero!(G)
    make_zero!(p.dspline_grid.eval)
    make_zero!(p.dspline_grid.basis_function_products)
    for spline_dimension in p.dspline_grid.spline_dimensions
        make_zero!(spline_dimension.eval)
    end
    autodiff(
        Reverse,
        image_loss,
        Active,
        Duplicated(control_points_flat, G),
        Duplicated(p.spline_grid, p.dspline_grid)
    )
    return nothing
end

x0 = zeros(prod(n_control_points))

dspline_grid = make_zero(spline_grid)
prob = OptimizationProblem(
    OptimizationFunction((x, p) -> image_loss(x, spline_grid); grad=image_loss_grad!),
    x0,
    (; spline_grid, dspline_grid),
)

# Optimize
sol = solve(prob, BFGS(); maxiters = 20)

# Plot solution
spline_grid.control_points .= reshape(sol.u, (n_control_points..., dim_out))
evaluate!(spline_grid)
plot(spline_grid)
```
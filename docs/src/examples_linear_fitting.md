# Linear fitting example

In this section we demonstrate how a spline grid can be fitted. We will fit a spline surface to the following image.

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

## Fitting

We fit the spline surface to the image by interpreting the spline grid evaluation as a linear mapping.

```@example tutorial
using LinearMaps
using IterativeSolvers
using Plots

spline_grid_map = LinearMap(spline_grid)
sol = lsqr(Matrix(spline_grid_map), vec(image_array)) # TODO: Stop requiring the matrix explicitly
spline_grid.control_points .= reshape(sol, size(spline_grid.control_points))
evaluate!(spline_grid)
plot(spline_grid)
```
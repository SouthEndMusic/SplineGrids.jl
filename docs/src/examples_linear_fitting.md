# Linear fitting example

In this section we demonstrate how a spline grid can be fitted. We will fit a spline surface to the following image.

```@example tutorial
using FileIO
image = load(normpath(@__DIR__, "julia_logo.png"))
```

## Defining the spline grid

We define a spline grid with 2 input dimensions, 1 output dimension and a sample grid which matches the image resolution.

```@example tutorial
using Colors: Gray
using SplineGrids

degree = (2, 2)
image_array = Float32.(Gray.(image[end:-1:1, :]))'
n_sample_points = size(image_array)
n_control_points = ntuple(i -> n_sample_points[i]÷10, 2)
dim_out = 1

spline_dimensions = SplineDimension.(n_control_points, degree, n_sample_points)
spline_grid = SplineGrid(spline_dimensions, dim_out)
spline_grid
```

The basis of this spline geometry looks like this:

```@example tutorial
using CairoMakie

plot_basis(spline_grid; title = "Unrefined spline basis")
```

## Fitting

We fit the spline surface to the image in a least squares sense, by interpreting the spline grid evaluation as a linear mapping.

```@example tutorial
using LinearMaps
using IterativeSolvers
using Plots

spline_grid_map = LinearMap(spline_grid)
sol = lsqr(spline_grid_map, vec(image_array))
copyto!(spline_grid.control_points, sol)
evaluate!(spline_grid)
Plots.plot(spline_grid; title = "Least squares fit")
```

## Matrix

The least-squares fitting procedure above is matrix free, but the linear mapping can be converted into a (sparse) matrix for inspection.

```@example tutorial
using SparseArrays

n_control_points = (5, 5)
degree = (2, 2)
n_sample_points = (15, 15)
dim_out = 1

spline_dimensions = SplineDimension.(n_control_points, degree, n_sample_points)
spline_grid_ = SplineGrid(spline_dimensions, dim_out)
spline_grid_map = LinearMap(spline_grid_)
M = sparse(spline_grid_map)
Plots.heatmap(M[end:-1:1,:])
```

## Local refinement informed by local error

Clearly the error of the fit is largest around the boundary of the text:

```@example tutorial
err_unrefined = (spline_grid.eval - image_array).^2
Plots.heatmap(err_unrefined[:,:,1]', colormap =  c=cgrad(:RdYlGn, rev=true))
title!("Squared error per pixel")
```

We can easily locally refine the spline basis by mapping this error back on to the control points.

```@example tutorial
spline_grid = add_default_local_refinement(spline_grid)
error_informed_local_refinement!(spline_grid, err_unrefined)
deactivate_overwritten_control_points!(spline_grid.control_points)
plot_basis(spline_grid; title = "Refined spline basis")
```

We can now fit the image again with the refined basis.

```@example tutorial
spline_grid_map = LinearMap(spline_grid)
sol = lsqr(spline_grid_map, vec(image_array))
copyto!(
    spline_grid.control_points, 
    reshape(sol, get_n_control_points(spline_grid.control_points), dim_out)
)
evaluate!(spline_grid.control_points)
evaluate!(spline_grid)
Plots.plot(spline_grid; title = "Least squares fit")
```

and the local error looks a lot better:

```@example tutorial
err_refined = (spline_grid.eval - image_array).^2
Plots.heatmap(err_refined[:,:,1]', colormap =  c=cgrad(:RdYlGn, rev=true), clims = (0, maximum(err_unrefined)))
title!("Squared error per pixel")
```
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

degree = (2, 2)
image_array = Float32.(Gray.(image[end:-1:1, :]))'
n_sample_points = size(image_array)
n_control_points = ntuple(i -> n_sample_points[i]÷10, 2)
dim_out = 1

spline_dimensions = SplineDimension.(n_control_points, degree, n_sample_points)
spline_grid = SplineGrid(spline_dimensions, dim_out)
spline_grid
```

## Fitting

We fit the spline surface to the image in a least squares sense, by interpreting the spline grid evaluation as a linear mapping.

```@example tutorial
using LinearMaps
using IterativeSolvers
using Plots

spline_grid_map = LinearMap(spline_grid)
sol = lsqr(spline_grid_map, vec(image_array))
spline_grid.control_points .= reshape(sol, size(spline_grid.control_points))
evaluate!(spline_grid)
plot(spline_grid)
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
heatmap(M[end:-1:1,:])
```

## Local refinement informed by local error

Clearly the error of the fit is largest around the boundary of the text:

```@example tutorial
err = (spline_grid.eval - image_array).^2
heatmap(err[:,:,1]', colormap =  c=cgrad(:RdYlGn, rev=true))
title!("Squared error per pixel")
```

We can easily locally refine the spline basis by mapping this error back on to the control points.

```@example tutorial
using CairoMakie

spline_grid = add_default_local_refinement(spline_grid)
error_informed_local_refinement!(spline_grid, err)
deactivate_overwritten_control_points!(spline_grid.control_points)
plot_basis(spline_grid)
```

We can now fit the image again with the refined basis (coming up).
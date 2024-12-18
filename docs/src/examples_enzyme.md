# Enzyme example

To demonstrate how Enzyme can be used to incorporate a spline grid into an optimization pipeline, we show what the derivative of a spline grid with respect to one control point looks like.

## Defining the spline grid

We define a spline grid with 2 input dimensions and 1 output dimension.

```@example tutorial
using SplineGrids

n_control_points = (10, 10)
degree = (2, 2)
n_sample_points = (100, 100)
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
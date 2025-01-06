# NURBS example

NURBS allow you to specify a weight for each control point, giving an extra degree of freedom. The spline grid is non-linear in terms of the weights. NURBS allow you to specify certain geometries exactly that splines can only approximate, like conic sections. The example below shows an exact representation of the unit circle.

```@example tutorial
using SplineGrids
using Plots

n_control_points = 9
degree = 2
n_sample_points = 500
dim_out = 2

knot_values = Float32[0, π / 2, π, 3π / 2, 2π]
multiplicities = Int32[3, 2, 2, 2, 3]
knot_vector = KnotVector(knot_values, multiplicities)

spline_dimension = SplineDimension(n_control_points, degree, n_sample_points; knot_vector)
nurbs_grid = NURBSGrid(spline_dimension, dim_out)
nurbs_grid.weights[2:2:end] .= 1 / sqrt(2)
nurbs_grid.control_points .= [
    1 0;
    1 1;
    0 1;
    -1 1;
    -1 0;
    -1 -1;
    0 -1;
    1 -1;
    1 0
]
evaluate!(nurbs_grid)
plot(nurbs_grid; aspect_ratio=:equal)
```
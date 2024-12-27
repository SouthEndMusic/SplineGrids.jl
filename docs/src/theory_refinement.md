# Refinement

The geometric operations discussed here are based on the following observation: When a knot is added to a knot vector, the basis functions defined on the original knot vector can be expressed as a linear combination of the basis functions defined on the new knot vector. This can be used to insert knots and with that new control points while the resulting spline geometry stays the same.

## Knot insertion

When a single knot is added to a knot vector this is called _knot insertion_. Let's have a look at what happens when we add a single knot.

```@example tutorial
using SplineGrids
using Random
using Plots

Random.seed!(2)

n_control_points = 8
degree = 2
n_sample_points = 500
Nout = 2

spline_dimension = SplineDimension(n_control_points, degree, n_sample_points; distribution=:random)
spline_grid = SplineGrid(spline_dimension, Nout)
for (i,θ) in enumerate(range(3π, 0, length = n_control_points))
    r = 2θ
    spline_grid.control_points[i, 1] = r * cos(θ)
    spline_grid.control_points[i, 2] = r * sin(θ)
end

plot(spline_dimension, title = "Original basis")
```

```@example tutorial
spline_grid_new = deepcopy(spline_grid)
spline_grid_new, refinement_matrix = insert_knot!(spline_grid_new, 0.25, 1)
spline_dimension_new = only(spline_grid_new.spline_dimensions)
plot(spline_dimension_new, title = "Basis after knot insertion")
```

Using the refinement matrix, we can for instance express the original basis function 3 in terms of the new basis functions 3 and 4.

```@example tutorial
data = decompress(spline_dimension)
data_new = decompress(spline_dimension_new)

(; sample_points) = spline_dimension
plot(sample_points, data[:, 3], label = "Original basis function 3")
plot!(sample_points, refinement_matrix[3, 3] * data_new[:, 3], label = "Scaled new basis function 3", ls = :dash)
plot!(sample_points, refinement_matrix[4, 3] * data_new[:, 4], label = "Scaled new basis function 4", ls = :dash)
ylims!(0,1)
```

If we define a curve in $\mathbb{R}^2$ with these bases, the refinement looks like this.

```@example tutorial
evaluate!(spline_grid)
plot(spline_grid, title = "Curve with original basis")
xlims!(-20, 15)
ylims!(-10, 20)
```

```@example tutorial
evaluate!(spline_grid_new)
plot(spline_grid_new, title = "Curve with basis after knot insertion")
xlims!(-20, 15)
ylims!(-10, 20)
```

## Knot refinement

We can also add multiple knots at a time, which is called _knot refinement_. A set of new knots can be supplied, but the default refinement behavior is to bisect each (non trivial) knot span.

```@example tutorial
spline_grid_new, refinement_matrix = refine!(spline_grid_new, 1)
evaluate!(spline_grid_new)
plot(spline_grid_new; title = "Curve with basis after knot refinement")
xlims!(-20, 15)
ylims!(-10, 20)
```
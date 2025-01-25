# Local refinement

As we have seen in the previous section ([Refinement](@ref)), The number of control points can be increased without changing the geometry by using refinement matrices. The key concepts of local refinement are as follows:

 1. Set up a base grid of control points
 2. Extend the control point grid using refinement matrices
 3. Overwrite certain control points in one of the intermediate extended grids

!!! note "Understanding THB splines"
    
    THB-splines can be a tricky construct to understand when you encounter them for the first time. This page does not aim to give a full theoretical background, that can be for instance found in [[1]](https://www.ag.jku.at/pubs/2016gjkmss.pdf). Playing around with THB-splines in the context of this package might however give you intuition for them more quickly than reading the theory.

## An example

Local refinement is probably best understood using an example. We use a surface here, because local refinement of curves is trivial.

We create a surface, and have a look at its basis functions.

```@example tutorial
using SplineGrids
using CairoMakie

n_control_points = (6, 6)
degree = (2, 2)
n_sample_points = (500, 500)
dim_out = 3

spline_dimensions = SplineDimension.(n_control_points, degree, n_sample_points)
spline_grid = SplineGrid(spline_dimensions, dim_out)
plot_basis(spline_grid)
```

This basis looks rather 'Cartesian', and knot insertion can only refine the basis in whole rows and columns. So what we do is set up refinement matrices for both dimensions, use some of the old basis functions and use some of the new ones.

```@example tutorial
# Set up `LocallyRefinedControlPoints` with a refinement matrix for both dimensions
spline_grid = add_default_local_refinement(spline_grid)
```

```@example tutorial
# Activate control point ranges (by default in the last extended control point grid).
# Note that some of these ranges overlap, this is fine and no duplicate control points are created
activate_local_control_point_range!(spline_grid, 1:4, 1:6)
activate_local_control_point_range!(spline_grid, 1:6, 1:2)
activate_local_control_point_range!(spline_grid, 9:10, 7:10)
deactivate_overwritten_control_points!(spline_grid)
plot_basis(spline_grid)
```

Note that some of the original basis functions are completely gone, which means that their contribution to the final geometry is completely overwritten. The function `deactivate_overwritten_control_points` weeds out the control points associated with these overwritten basis functions. This means that every active control point is guaranteed to influence the spline geometry (assuming there is at least one global sample point in the effective support of the basis function associated with that control point).

A nice property of this construction is that it can be iterated, creating a hierarchy. Let's refine the basis some more:

```@example tutorial
spline_grid = add_default_local_refinement(spline_grid)
```

```@example tutorial
activate_local_control_point_range!(spline_grid, 5:12, 1:4)
activate_local_control_point_range!(spline_grid, 7:8, 5:6)
deactivate_overwritten_control_points!(spline_grid)
plot_basis(spline_grid)
```

This is in fact an exact reproduction of the THB-spline example shown in [[1]](https://www.ag.jku.at/pubs/2016gjkmss.pdf) (p. 6).

Finally, let's give the $z$-coordinates of the control points some random values and look at the surface.

```@example tutorial
using Plots
using Random
Random.seed!(42)

spline_grid.control_points[:, 3] .= rand(Float32, get_n_control_points(spline_grid))
evaluate!(spline_grid.control_points)
evaluate!(spline_grid)
p = Plots.plot(spline_grid, camera = (30, 60, 1.5))
Plots.zlims!(p, 0, 1)
```

See [Local refinement informed by local error](@ref) for an example of how a local fitting error can be used to inform where to refine.

## References

[1] Giannelli, C., Jüttler, B., Kleiss, S. K., Mantzaflaris, A., Simeon, B., & Špeh, J. (2016). _THB-splines: An effective mathematical technology for adaptive refinement in geometric design and isogeometric analysis._ Computer Methods in Applied Mechanics and Engineering, 299, 337-365.

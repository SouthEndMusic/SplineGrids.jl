# Linear fitting example

In this section we demonstrate how a spline grid can be fitted. We will fit a spline surface to the following image.

```@example tutorial
using FileIO
using CairoMakie
image = rotr90(load(normpath(@__DIR__, "julia_logo.png")))
fig = Figure()
ratio = size(image, 1) / size(image, 2)
ax = Axis(fig[1, 1], aspect = ratio)
CairoMakie.image!(ax, image)
fig
```

## Defining the spline grid

We define a spline grid with 2 input dimensions, 1 output dimension and a sample grid which matches the image resolution.

```@example tutorial
using Colors: Gray
using SplineGrids

degree = (2, 2)
image_array = Float32.(Gray.(image))
n_sample_points = size(image_array)
n_control_points = ntuple(i -> n_sample_points[i]÷25, 2)
extent = ntuple(i -> (0, n_sample_points[i]), 2)
dim_out = 1

spline_dimensions = ntuple(i -> SplineDimension(n_control_points[i], degree[i], n_sample_points[i]; extent = extent[i]), 2)
spline_grid = SplineGrid(spline_dimensions, dim_out)
spline_grid
```

The basis of this spline geometry looks like this:

```@example tutorial
fig_total = Figure(size = (1600, 900)) # hide
plot_basis!(fig_total, spline_grid; i = 1, j = 1, title = "Unrefined spline basis") # hide
plot_basis(spline_grid; title = "Unrefined spline basis")
```

## Fitting

We fit the spline surface to the image in a least squares sense, by interpreting the spline grid evaluation as a linear mapping.

```@example tutorial
using LinearMaps
using IterativeSolvers
using Plots

function fit!(spline_grid)
    spline_grid_map = LinearMap(spline_grid)
    sol = lsqr(spline_grid_map, vec(image_array))
    copyto!(
        spline_grid.control_points,
        reshape(sol, get_n_control_points(spline_grid.control_points), dim_out)
    )
    evaluate!(spline_grid.control_points)
    evaluate!(spline_grid)
end

function plot_fit(spline_grid)
    Plots.plot(spline_grid; aspect_ratio = :equal, title = "Least squares fit", clims = (-0.5, 1.5), cmap = :viridis)
end
# hide
function _plot_fit(spline_grid, j) # hide
    ax_fit = Axis(fig_total[2, j], aspect = ratio; title = "Least squares fit") # hide
    CairoMakie.heatmap!(ax_fit, spline_grid.eval[:, :, 1], colorrange = (-0.5, 1.5)) # hide
end # hide

fit!(spline_grid)
_plot_fit(spline_grid, 1) # hide
plot_fit(spline_grid)
```

## Matrix

The least-squares fitting procedure above is matrix free, but the linear mapping can be converted into a (sparse) matrix for inspection.

```@example tutorial
using SparseArrays

# Make an analogous spline geometry with a smaller sample grid so the output dimensionality is not too large
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

Clearly the error of the fit is largest where the text is:

```@example tutorial
err_unrefined = (spline_grid.eval - image_array).^2

function plot_error(error)
    Plots.heatmap(error[:,:,1]', colormap = cgrad(:RdYlGn, rev=true), aspect_ratio = :equal, clims = (0, 1))
    title!("Squared error per pixel")
end
# hide
function _plot_error(error, j) # hide
    ax_err = Axis(fig_total[3, j], aspect = ratio; title = "Squared error per pixel") # hide
    CairoMakie.heatmap!(ax_err, error[:, :, 1], colormap = cgrad(:RdYlGn, rev=true); colorrange = (0, 1)) # hide
end # hide

_plot_error(err_unrefined, 1) # hide
plot_error(err_unrefined)
```

We can easily locally refine the spline basis by mapping this error back on to the control points.

```@example tutorial
function refine(spline_grid)
    spline_grid = add_default_local_refinement(spline_grid)
    error_informed_local_refinement!(spline_grid, err_unrefined)
    deactivate_overwritten_control_points!(spline_grid.control_points)
    spline_grid
end

spline_grid = refine(spline_grid)
plot_basis!(fig_total, spline_grid; i = 1, j = 2, title = "Refined spline basis (level 1)") # hide
plot_basis(spline_grid; title = "Refined spline basis (level 1)")
```

We can now fit the image again with the refined basis.

```@example tutorial
fit!(spline_grid)
_plot_fit(spline_grid, 2) # hide
plot_fit(spline_grid)
```

and the local error looks a bit better:

```@example tutorial
err_refined = (spline_grid.eval - image_array).^2
_plot_error(err_refined, 2) # hide
plot_error(err_refined)
```

## Iterating local refinement

Let's iterate the local refinement and fitting procedure a few more times to get a nicer result!

```@example tutorial
function _iteration(spline_grid, j) # hide
    spline_grid = refine(spline_grid) # hide
    plot_basis!(fig_total, spline_grid; i = 1, j,  title = "Refined spline basis (level $(j - 1))") # hide
    fit!(spline_grid) # hide
    _plot_fit(spline_grid, j) # hide
    err_refined = (spline_grid.eval - image_array).^2 # hide
    _plot_error(err_refined, j) # hide
    spline_grid # hide
end # hide
spline_grid = _iteration(spline_grid, 3) # hide
_iteration(spline_grid, 4) # hide
fig_total # hide
```
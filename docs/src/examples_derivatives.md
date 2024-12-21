# Derivatives example

The code below shows an example of how derivatives can be used to compute the normal vectors of a surface.

### Defining the surface

```@setup tutorial
using Base.Threads

# Set the number of threads to the number of available CPU cores
ENV["JULIA_NUM_THREADS"] = string(nthreads())
```

```@example tutorial
using SplineGrids

# Spline Grid Parameters
n_control_points = (5, 5)
degree = (2, 2)
n_sample_points = (100, 100)
Nout = 3
max_derivative_order = 1

# Define the spline grid
spline_dimensions = SplineDimension.(n_control_points, degree, n_sample_points; max_derivative_order)
spline_grid = SplineGrid(spline_dimensions, Nout)
for i in 1:n_control_points[1]
    for j in 1:n_control_points[2]
        spline_grid.control_points[i, j, 3] = exp(-((i - 3)^2 + (j - 3)^2))
    end
end
```

### Inspecting the basis functions

Let's look at the basis functions for the first dimension:

```@example tutorial
using Plots

plot(spline_dimensions[1])
```

and the derivatives of these basis functions:

```@example tutorial
using Plots

plot(spline_dimensions[1]; derivative_order = 1)
```

### Evaluating the surface

The data is copied because `spline_grid.eval` will be overwritten below when computing the partial derivatives.

```@example tutorial
evaluate!(spline_grid)
spline_grid_data = copy(spline_grid.eval)
nothing
```

### Computing the partial derivatives

```@example tutorial
# Derivatives of the surface with respect to the first dimension parameter
evaluate!(spline_grid; derivative_order=(1, 0))
∂₁spline_grid_data = copy(spline_grid.eval)

# Derivatives of the surface with respect to the second dimension parameter
evaluate!(spline_grid; derivative_order=(0, 1))
∂₂spline_grid_data = spline_grid.eval
nothing
```

### Computing the normal vectors

```@example tutorial
using LinearAlgebra

normal_vectors = similar(spline_grid_data)
for i in 1:n_sample_points[1]
    for j in 1:n_sample_points[2]
        normal_vectors[i, j, :] .= cross(
            view(∂₁spline_grid_data, i, j, :),
            view(∂₂spline_grid_data, i, j, :)
        )
        normal_vectors[i, j, :] ./= 10 * norm(view(normal_vectors, i, j, :))
    end
end
```

### Plotting

```@example tutorial
using CairoMakie

# Plotting the surface
f = Figure()
ax = Axis3(f[1, 1])
CairoMakie.surface!(ax, eachslice(spline_grid_data, dims=3)...)

# Plotting a subset of the normal vectors
spline_grid_data_reshaped = reshape(spline_grid_data, (prod(n_sample_points), 3))
normal_vectors_reshaped = reshape(normal_vectors, (prod(n_sample_points), 3))

CairoMakie.quiver!(
    ax,
    eachslice(view(spline_grid_data_reshaped, 1:192:prod(n_sample_points), :), dims=2)...,
    eachslice(view(normal_vectors_reshaped, 1:192:prod(n_sample_points), :), dims=2)...,
    arrowsize=0.03,
)

CairoMakie.xlims!(0, 1)
CairoMakie.ylims!(0, 1)

f
```
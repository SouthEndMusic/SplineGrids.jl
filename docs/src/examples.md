# Examples

This page shows examples of spline grids of various dimensionalities via their plotting recipes.

## 1 input, 1 output

```@example tutorial
using SplineGrids
using Plots

# Spline grid parameters
n_control_points = 10
degree = 3
n_sample_points = 250
dim_out = 1

# Define spline grid
spline_dimension = SplineDimension(n_control_points, degree, n_sample_points, extent=(3.0, 5.0))
spline_grid = SplineGrid(spline_dimension, dim_out)

# Set control points
spline_grid.control_points .= [0.342, 0.633, 0.446, 0.716, 0.843, 0.171, 0.061, 0.973, 0.057, 0.671]

# Evaluate spline grid
evaluate!(spline_grid)

# Plot
plot(spline_grid)
```

```@example tutorial
spline_grid
```

## 1 input, 2 outputs

```@example tutorial
# Spline grid parameters
n_control_points = 10
degree = 2
n_sample_points = 250
dim_out = 2

# Define spline grid
spline_dimension = SplineDimension(n_control_points, degree, n_sample_points)
spline_grid = SplineGrid(spline_dimension, dim_out)

# Set control points
r = 1:n_control_points
θ = 1:n_control_points
spline_grid.control_points[:, 1] .= @. r * cos(θ)
spline_grid.control_points[:, 2] .= @. r * sin(θ)

# Evaluate spline grid
evaluate!(spline_grid)

# Plot
plot(spline_grid)
```

```@example tutorial
spline_grid
```

## 1 input, 3 outputs

```@example tutorial
# Spline grid parameters
n_control_points = 10
degree = 2
n_sample_points = 250
dim_out = 3

# Define spline grid
spline_dimension = SplineDimension(n_control_points, degree, n_sample_points)
spline_grid = SplineGrid(spline_dimension, dim_out)

# Set control points
r = 1:n_control_points
θ = 1:n_control_points
z = 1:n_control_points
spline_grid.control_points[:, 1] .= @. r * cos(θ)
spline_grid.control_points[:, 2] .= @. r * sin(θ)
spline_grid.control_points[:, 3] .= z

# Evaluate spline grid
evaluate!(spline_grid)

# Plot
plot(spline_grid)
```

```@example tutorial
spline_grid
```

## 1 input, 4 outputs

```@example tutorial
# Spline grid parameters
n_control_points = 10
degree = 2
n_sample_points = 250
dim_out = 4

# Define spline grid
spline_dimension = SplineDimension(n_control_points, degree, n_sample_points)
spline_grid = SplineGrid(spline_dimension, dim_out)

# Set control points
r = 1:n_control_points
θ = 1:n_control_points
z = 1:n_control_points
c = sin.(z)
spline_grid.control_points[:, 1] .= @. r * cos(θ)
spline_grid.control_points[:, 2] .= @. r * sin(θ)
spline_grid.control_points[:, 3] .= z
spline_grid.control_points[:, 4] .= c

# Evaluate spline grid
evaluate!(spline_grid)

# Plot
plot(spline_grid)
```

```@example tutorial
spline_grid
```

## 2 inputs, 1 output

```@example tutorial
# Spline grid parameters
n_control_points = (5, 6)
degree = (2, 2)
n_sample_points = (60, 50)
dim_out = 1

# Define spline grid
spline_dimensions = SplineDimension.(n_control_points, degree, n_sample_points)
spline_grid = SplineGrid(spline_dimensions, dim_out)

# Set control points
spline_grid.control_points .= 0
spline_grid.control_points[1:2:end] .= 1:prod(n_control_points)/2

# Evaluate spline grid 
evaluate!(spline_grid)

# Plot
plot(spline_grid)
```

```@example tutorial
spline_grid
```

## 2 inputs, 2 outputs

```@example tutorial
# Spline grid parameters
n_control_points = (4, 4)
degree = (2, 3)
n_sample_points = (15, 20)
dim_out = 2

# Define spline grid
spline_dimensions = SplineDimension.(n_control_points, degree, n_sample_points)
spline_grid = SplineGrid(spline_dimensions, dim_out)

# Set control points
spline_grid.control_points .= [-0.358 0.795 -0.016 -0.295; 0.923 -0.182 -0.644 0.612; -0.91 0.708 -0.426 0.412; 0.56 -0.78 0.515 0.676;;; 0.938 0.393 -0.702 -0.99; -0.578 0.305 -0.842 -0.57; 0.034 -0.813 -0.514 0.162; -0.016 -0.822 -0.261 -0.148]

# Evaluate spline grid
evaluate!(spline_grid)

# plot
plot(spline_grid)
```

```@example tutorial
spline_grid
```

## 2 inputs, 3 outputs

```@example tutorial
# Spline grid parameters
n_control_points = (25, 2)
degree = (3, 1)
n_sample_points = (100, 100)
dim_out = 3

# Define spline grid
spline_dimension = SplineDimension.(n_control_points, degree, n_sample_points)
spline_grid = SplineGrid(spline_dimension, dim_out)

# Set control points
spline_grid.control_points .= 0
R = 3
r = 1
ρ = range(-r, r, length=n_control_points[2])
for (i, θ) in enumerate(range(0, 2π, length=n_control_points[1]))
    ϕ = 2θ
    spline_grid.control_points[i, :, 1] .= @. (R + ρ * cos(ϕ)) * cos(θ)
    spline_grid.control_points[i, :, 2] .= @. (R + ρ * cos(ϕ)) * sin(θ)
    spline_grid.control_points[i, :, 3] .= ρ * sin(ϕ)
end

# Evaluate spline grid
evaluate!(spline_grid)

# Plot
plot(spline_grid)
```

```@example tutorial
spline_grid
```
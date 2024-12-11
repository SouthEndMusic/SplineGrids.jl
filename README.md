# SplineGrids.jl

[![Global Docs](https://img.shields.io/badge/docs-blue.svg)](https://southendmusic.github.io/SplineGrids.jl/)
[![CI](https://github.com/SouthEndMusic/SplineGrids.jl/actions/workflows/Tests.yml/badge.svg?branch=master)](https://github.com/SouthEndMusic/SplineGrids.jl/actions/workflows/Tests.yml)

SplineGrids.jl is designed to efficiently evaluate a broad class of spline objects on a grid in the spline's domain. The package supports:
- Any number of input and output dimensions (see the examples [here](https://southendmusic.github.io/SplineGrids.jl/dev/examples_dimensions/))
- Any degree of basis functions and type of knot vector
- Any combination of partial derivatives
- Using weights to define NURBS (coming up)
- Local refinement (coming up)
- CPU and GPU backends via [KernelAbstractions.jl](https://github.com/JuliaGPU/KernelAbstractions.jl)

The intended use case of this package is to be integrated into the [SciML ecosystem](https://sciml.ai/), for easy and efficient incorporation of spline geometries in problems of fitting, optimization, differential equation solving and machine learning. Since SplineGrids.jl makes heavy use of in-place computations, the recommended automatic differentiation framework to use with SplineGrids.jl is [Enzyme](https://github.com/EnzymeAD/Enzyme.jl).

## API

```julia
using SplineGrids

# Parameters per input dimension
n_control_points = (10, 10, 5)
degree = (2, 3, 2)
n_sample_points = (50, 50, 25)

# The number of output dimensions
Nout = 4

# Defining the spline grid
spline_dimensions = SplineDimension.(n_control_points, degree, n_sample_points)
spline_grid = SplineGrid(spline_dimensions, Nout)

# Set the desired control points
spline_grid.control_points .= rand(n_control_points..., Nout)

# Evaluate
evaluate!(spline_grid)

# The output can be found in spline_grid.eval of shape (n_sample_points..., Nout)
```

## History and theory

The most well known introduction to spline theory is probably _The NURBS book_ [1].

SplineGrids.jl was inspired by work on caustic design by optimizing spline-based lens surfaces with differentiable ray tracing. This work was published in a [master thesis](https://resolver.tudelft.nl/uuid:0c514716-f2db-455e-b75d-3cf9cfeed8bb) and later in follow-up research [2][3][4].

Some of the core ideas for this package where implemented earlier in [NURBS Pytorch](https://github.com/SouthEndMusic/NURBS_PyTorch), but that package was never properly tested or released.

## References

[1] Piegl, L., & Tiller, W. (2012). _The NURBS book_. Springer Science & Business Media.

[2] Koning, B. D., Heemels, A., Adam, A., & Möller, M. (2024). _Gradient descent-based freeform optics design for illumination using algorithmic differentiable non-sequential ray tracing_. Optimization and Engineering, 25(3), 1203-1235.

[3] Heemels, A., De Koning, B., Möller, M., & Adam, A. (2024). _Optimizing freeform lenses for extended sources with algorithmic differentiable ray tracing and truncated hierarchical B-splines_. Optics Express, 32(6), 9730-9746.

[4] Heemels, A., de Koning, B., Moller, M., & Adama, A. (2024, May). _Unsupervised design of illumination optics using algorithmic differentiable raytracing_. In International Optical Design Conference 2023 (Vol. 12798, pp. 84-85). SPIE.

abstract type AbstractSplineGrid{Nin, Nout} end

"""
The SplineGrid is the central object of the `SplineGrids.jl` package, containing
all information to evaluate the defined spline on the defined grid.

## Fields

  - `spline_dimensions`: A SplineDimension per dimension of the spline, containing data to evaluate
    basis functions.
  - `control points`: The points that define the shape of the spline, and in how many dimensions it is embedded.
  - `weights`: For now unsupported, will eventually be used to define NURBS.
  - `eval`: The array where the evaluated spline is stored.
  - `basis_function_products`: An array of intermediate results for evaluating the spline grids, containing products of basis functions
    from the various spline dimensions.
"""
struct SplineGrid{
    S <: SplineDimension, C <: AbstractArray, W <: Union{AbstractArray, Nothing},
    E <: AbstractArray, B <: AbstractArray, Nin, Nout} <: AbstractSplineGrid{Nin, Nout}
    spline_dimensions::NTuple{Nin, S}
    control_points::C
    weights::W
    eval::E
    basis_function_products::B
    function SplineGrid(
            spline_dimensions, control_points, weights, eval, basis_function_products)
        # TODO: Add validation of combination of control points, weights, and basis functions
        Nin = length(spline_dimensions)
        Nout = size(control_points)[end]
        new{eltype(spline_dimensions), typeof(control_points), typeof(weights),
            typeof(eval), typeof(basis_function_products), Nin, Nout}(
            spline_dimensions, control_points, weights, eval, basis_function_products)
    end
end

# TODO: Add indexing functionality such that indexing a spline_grid implies indexing spline_grid.eval

"""
    SplineGrid(spline_dimensions::NTuple{N_in, <:SplineDimension}, dim_out::Integer)

Define a `SplineGrid` from an NTuple of spline dimensions and the number of output dimensions.

## Inputs

  - `spline_dimensions`: an NTuple of spline dimensions
  - `dim_out`: The number of output dimensions. I.e. the control points and thus the spline live in ℝ^dim_out.
"""
function SplineGrid(spline_dimensions::NTuple{N_in, <:SplineDimension},
        dim_out::Integer)::SplineGrid where {N_in}
    # The size of the point grid on which the spline is evaluated
    size_eval_grid = (length(spline_dimension.sample_points) for spline_dimension in spline_dimensions)
    # The size of the grid of control points
    size_cp_grid = get_n_basis_functions.(spline_dimensions)
    # The actual control points
    control_points = zeros(size_cp_grid..., dim_out)
    set_unit_cp_grid!(control_points)
    # Preallocated memory for basis function product evaluation
    basis_function_products = zeros(size_eval_grid...)
    # Preallocated memory for grid evaluation of the spline
    eval = zeros(size_eval_grid..., dim_out)
    # NURBS are not supported yet
    weights = nothing
    SplineGrid(spline_dimensions, control_points, weights, eval, basis_function_products)
end

# Support inputting a single spline dimension instead of a tuple
function SplineGrid(spline_dimension::SplineDimension, args...; kwargs...)
    SplineGrid((spline_dimension,), args...; kwargs...)
end

# Get the size of the block of control points that each output of the spline
# depends on
function cp_kernel_size(spline_grid::SplineGrid)
    Tuple(spline_dim.degree + 1 for spline_dim in spline_grid.spline_dimensions)
end

# Outer product of n vectors, with thanks to Michael Abbott
function outer!(A::AbstractArray{T, N}, vs::Vararg{SubArray, N}) where {T, N}
    vecs = ntuple(n -> reshape(vs[n], ntuple(Returns(1), n - 1)..., :), N)
    broadcast!(*, A, vecs...)
end

"""
    evaluate!(spline_grid::SplineGrid)

Evaluate the spline grid, that is: take the evaluated basis functions for each sample point
for each SplineDimension, and compute the output grid on each sample point combination
as a linear combination of control with basis function products as coefficients.
"""
function evaluate!(spline_grid::SplineGrid)::Nothing
    (; basis_function_products, eval, spline_dimensions, control_points) = spline_grid
    eval .= 0

    # Loop over the positions in the kernel of control points each spline evaluation depends on
    for I in CartesianIndices(cp_kernel_size(spline_grid))
        # Compute basis function products as an outer product of the basis function values per 
        # spline dimension
        outer!(basis_function_products,
            (view(spline_dim.eval, :, i) for (i, spline_dim) in zip(
                Tuple(I), spline_dimensions))...)
        # Loop over all sample points to compute one one basis_function_product * control_point contribution
        for J in CartesianIndices(basis_function_products)
            cp_indices = (spline_dim.sample_indices[j] - spline_dim.degree - 1 + i for (i, j, spline_dim) in zip(
                Tuple(I), Tuple(J), spline_dimensions))
            control_point = view(control_points, cp_indices..., :)
            @. eval[J, :] += basis_function_products[J] * control_point
        end
    end
    return nothing
end
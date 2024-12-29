"""
The SplineGrid is the central object of the `SplineGrids.jl` package, containing
all information to evaluate the defined spline on the defined grid.

## Fields

  - `spline_dimensions`: A SplineDimension per dimension of the spline, containing data to evaluate
    basis functions.
  - `control points`: The points that define the shape of the spline, and in how many dimensions it is embedded.
  - `denominator`: An optional array of intermediate results to evaluate the denominator of NURBS.
  - `weights`: Control point weights to define NURBS.
  - `eval`: The array where the evaluated spline grid is stored.
  - `sample_indices`: For each global sample point, the linear index in the `control_points` array before it is offset
    for a particular index in the control point kernel and output dimension.
  - `basis_function_products`: An array of intermediate results for evaluating the spline grids, containing products of basis functions
    from the various spline dimensions.
"""
struct SplineGrid{
    S <: AbstractSplineDimension{T} where {T},
    C <: AbstractArray{T} where {T}, # N = Nin + 1
    D <: Union{AbstractArray{T, Nin}, Nothing} where {T, Nin},
    W <: Union{AbstractArray{T, Nin}, Nothing} where {T, Nin},
    E <: AbstractArray{T} where {T},
    I <: AbstractArray{<:Integer, Nin} where {Nin},
    B <: AbstractArray{T} where {T},
    Nin,
    Nout,
    HasWeights,
    T <: AbstractFloat
} <: AbstractSplineGrid{Nin, Nout, HasWeights, T}
    spline_dimensions::NTuple{Nin, S}
    control_points::C
    denominator::D
    weights::W
    eval::E
    sample_indices::I
    basis_function_products::B
    function SplineGrid(
            spline_dimensions,
            control_points,
            eval,
            sample_indices,
            basis_function_products;
            denominator = nothing,
            weights = nothing)
        # TODO: Add validation of combination of control points, weights, and basis functions
        new{
            eltype(spline_dimensions),
            typeof(control_points),
            typeof(denominator),
            typeof(weights),
            typeof(eval),
            typeof(sample_indices),
            typeof(basis_function_products),
            length(spline_dimensions),
            size(control_points)[end],
            isa(weights, AbstractArray),
            eltype(control_points)}(
            spline_dimensions,
            control_points,
            denominator,
            weights,
            eval,
            sample_indices,
            basis_function_products
        )
    end
end

# TODO: Add indexing functionality such that indexing a spline_grid implies indexing spline_grid.eval

"""
    SplineGrid(spline_dimensions::NTuple{Nin, <:SplineDimension}, Nout::Integer)

Define a `SplineGrid` from an NTuple of spline dimensions and the number of output dimensions.

## Inputs

  - `spline_dimensions`: an NTuple of spline dimensions
  - `Nout`: The number of output dimensions. I.e. the control points and thus the spline live in ℝ^Nout.
"""
function SplineGrid(
        spline_dimensions::NTuple{Nin, <:AbstractSplineDimension{T}},
        Nout::Integer
)::AbstractSplineGrid{Nin} where {Nin, T}
    backend = get_backend(first(spline_dimensions))
    # TODO: Then do the same for NURBSGrid
    # TODO: Then look for other places where arrays are allocated with zeros/ones 

    # The size of the point grid on which the spline is evaluated
    size_eval_grid = get_sample_grid_size(spline_dimensions)
    # The size of the grid of control points
    size_cp_grid = get_n_basis_functions.(spline_dimensions)
    # The control points
    control_points = allocate(backend, T, size_cp_grid..., Nout)
    set_unit_cp_grid!(control_points)
    # Preallocated memory for basis function product evaluation
    basis_function_products = allocate(backend, T, size_eval_grid...)
    # Preallocated memory for grid evaluation of the spline
    eval = allocate(backend, T, size_eval_grid..., Nout)
    # Linear indices for control points per global sample point
    sample_indices = get_global_sample_indices(spline_dimensions, control_points)
    SplineGrid(
        spline_dimensions,
        control_points,
        eval,
        sample_indices,
        basis_function_products
    )
end

# Support inputting a single spline dimension instead of a tuple
function SplineGrid(spline_dimension::SplineDimension, args...; kwargs...)
    SplineGrid((spline_dimension,), args...; kwargs...)
end

@kernel function spline_muladd_kernel(
        eval,
        @Const(basis_function_products),
        @Const(control_points),
        @Const(sample_indices),
        offset
)
    # Index of the global sample point
    J = @index(Global, Cartesian)

    # Output dimensionality
    Nout = size(control_points)[end]

    # The total number of control points
    cp_grid_size = prod(size(control_points)[1:(end - 1)]) # Could be done outside kernel

    # The product of basis functions for the current sample point
    # and the current control point kernel location
    basis_function_product = basis_function_products[J]

    # The linear index of the required control point
    lin_cp_index_base = sample_indices[J] + offset

    # The contribution of the required control point
    for dim_out in 1:Nout
        lin_cp_idx = lin_cp_index_base + dim_out * cp_grid_size
        eval[J, dim_out] += basis_function_product * control_points[lin_cp_idx]
    end
end

"""
    evaluate!(spline_grid::AbstractSplineGrid{Nin};
        derivative_order::NTuple{Nin, <:Integer} = ntuple(_ -> 0, Nin),
        control_points::AbstractArray = spline_grid.control_points,
        eval::AbstractArray = spline_grid.eval)

Evaluate the spline grid, that is: take the evaluated basis functions for each sample point
for each SplineDimension, and compute the output grid on each sample point combination
as a linear combination of control with basis function products as coefficients.

Uses the `control_points` and `eval` arrays from the `spline_grid` by default,
but different arrays can be specified as a convenience for optimization algorithms.
"""
function evaluate!(spline_grid::AbstractSplineGrid{Nin};
        derivative_order::NTuple{Nin, <:Integer} = ntuple(_ -> 0, Nin),
        control_points::AbstractArray = spline_grid.control_points,
        eval::AbstractArray = spline_grid.eval
)::Nothing where {Nin}
    (; basis_function_products, spline_dimensions, sample_indices) = spline_grid
    validate_partial_derivatives(spline_dimensions, derivative_order)
    eval .= 0
    @assert size(control_points) == size(spline_grid.control_points)
    @assert size(eval) == size(spline_grid.eval)

    control_point_kernel_size = get_cp_kernel_size(spline_dimensions)
    backend = get_backend(eval)
    kernel! = spline_muladd_kernel(backend)

    # Loop over the positions in the kernel of control points each spline evaluation depends on
    for I in CartesianIndices(control_point_kernel_size)
        # Compute basis function products as an outer product of the basis function values per 
        # spline dimension
        outer!(
            basis_function_products,
            (
                view(spline_dim.eval, :, i, derivative_order_ + 1) for
            (i, spline_dim, derivative_order_) in zip(
                Tuple(I), spline_dimensions, derivative_order))...
        )

        # Add the 'basis function product * control point' contribution to eval
        offset = get_offset(size(control_points), Tuple(I))
        kernel!(
            eval,
            basis_function_products,
            control_points,
            sample_indices,
            offset,
            ndrange = size(sample_indices))
        synchronize(backend)
    end
    return nothing
end

@kernel function spline_muladd_adjoint_kernel(
        control_points,
        @Const(eval),
        @Const(basis_function_products),
        @Const(sample_indices),
        offset
)
    # Index of the global sample point
    J = @index(Global, Cartesian)

    # Output dimensionality
    Nout = size(control_points)[end]

    # The total number of control points
    cp_grid_size = prod(size(control_points)[1:(end - 1)]) # Could be done outside kernel

    # The product of basis functions for the current sample point
    # and the current control point kernel location
    basis_function_product = basis_function_products[J]

    # The linear index of the required control point
    lin_cp_index_base = sample_indices[J] + offset

    for dim_out in 1:Nout
        lin_cp_idx = lin_cp_index_base + dim_out * cp_grid_size
        Atomix.@atomic control_points[lin_cp_idx] += basis_function_product *
                                                     eval[J, dim_out]
    end
end

"""
    evaluate_adjoint!(spline_grid::AbstractSplineGrid{Nin};
        derivative_order::NTuple{Nin, <:Integer} = ntuple(_ -> 0, Nin),
        control_points::AbstractArray = spline_grid.control_points,
        eval::AbstractArray = spline_grid.eval)

evaluate the adjoint of the linear mapping `control_points -> eval`. This is a computation of the form
`eval -> control_points`. If we write `evaluate!(spline_grid)` as a matrix vector multiplication `eval = M * control_points`,
Then the adjoint is given by `v -> M' * v`. This mapping is used in fitting algorithms.
"""
function evaluate_adjoint!(spline_grid::AbstractSplineGrid{Nin};
        derivative_order::NTuple{Nin, <:Integer} = ntuple(_ -> 0, Nin),
        control_points::AbstractArray = spline_grid.control_points,
        eval::AbstractArray = spline_grid.eval
)::Nothing where {Nin}
    (; basis_function_products, spline_dimensions, sample_indices) = spline_grid
    validate_partial_derivatives(spline_dimensions, derivative_order)
    control_points .= 0
    @assert size(control_points) == size(spline_grid.control_points)
    @assert size(eval) == size(spline_grid.eval)

    control_point_kernel_size = get_cp_kernel_size(spline_dimensions)
    backend = get_backend(eval)
    kernel! = spline_muladd_adjoint_kernel(backend)

    # Loop over the positions in the kernel of control points each spline evaluation depends on
    for I in CartesianIndices(control_point_kernel_size)
        # Compute basis function products as an outer product of the basis function values per 
        # spline dimension
        outer!(
            basis_function_products,
            (
                view(spline_dim.eval, :, i, derivative_order_ + 1) for
            (i, spline_dim, derivative_order_) in zip(
                Tuple(I), spline_dimensions, derivative_order))...
        )

        # Add the 'basis function product * control point' contribution to eval
        offset = get_offset(size(control_points), Tuple(I))
        kernel!(control_points, eval, basis_function_products,
            sample_indices, offset, ndrange = size(sample_indices))
        synchronize(backend)
    end
    return nothing
end
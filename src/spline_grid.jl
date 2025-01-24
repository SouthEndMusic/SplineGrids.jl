"""
    SplineGrid(
        spline_dimensions,
        control_points,
        weights,
        eval)

The SplineGrid is the central object of the `SplineGrids.jl` package, containing
all information to evaluate the defined spline on the defined grid.

## Fields

  - `spline_dimensions`: A SplineDimension per dimension of the spline, containing data to evaluate
    basis functions.
  - `control points`: The points that define the shape of the spline, and in how many dimensions it is embedded.
  - `weights`: Control point weights to define NURBS.
  - `eval`: The array where the evaluated spline grid is stored.
"""
struct SplineGrid{
    Nin,
    Nout,
    Tv <: AbstractFloat,
    Ti <: Integer,
    S <: SplineDimension{Tv, Ti},
    C <: AbstractControlPoints{Nin, Nout, Tv},
    W <: Union{AbstractArray{Tv, Nin}, Nothing},
    E <: AbstractArray{Tv},
    HasWeights
} <: AbstractSplineGrid{Nin, Nout, HasWeights, Tv, Ti}
    spline_dimensions::NTuple{Nin, S}
    control_points::C
    weights::W
    eval::E
    function SplineGrid(
            spline_dimensions,
            control_points,
            weights,
            eval
    )
        validate_spline_grid(spline_dimensions, control_points, weights, eval)
        new{
            length(spline_dimensions),
            size(control_points)[end],
            eltype(control_points),
            eltype(first(spline_dimensions).sample_indices),
            eltype(spline_dimensions),
            typeof(control_points),
            typeof(weights),
            typeof(eval),
            isa(weights, AbstractArray)
        }(
            spline_dimensions,
            control_points,
            weights,
            eval
        )
    end
end

"""
    SplineGrid(spline_dimensions::NTuple{Nin, <:SplineDimension{Tv, Ti}}, Nout::Integer)::SplineGrid{Nin, Tv, Ti} where {Nin, Tv, Ti}

Define a `SplineGrid` from an NTuple of spline dimensions and the number of output dimensions.

## Inputs

  - `spline_dimensions`: an NTuple of spline dimensions
  - `Nout`: The number of output dimensions. I.e. the control points and thus the spline live in ℝ^Nout.
"""
function SplineGrid(
        spline_dimensions::NTuple{Nin, <:SplineDimension{Tv, Ti}},
        Nout::Integer
)::SplineGrid{Nin, Nout, Tv, Ti} where {Nin, Tv, Ti}
    backend = get_backend(first(spline_dimensions))
    # The size of the point grid on which the spline is evaluated
    size_eval_grid = get_sample_grid_size(spline_dimensions)
    # The size of the grid of control points
    size_cp_grid = get_control_point_grid_size(spline_dimensions)
    # The control points
    control_points = KernelAbstractions.zeros(CPU(), Tv, size_cp_grid..., Nout)
    set_unit_cp_grid!(control_points)
    control_points = DefaultControlPoints(control_points)
    control_points = adapt(backend, control_points)
    # Preallocated memory for grid evaluation of the spline
    eval = KernelAbstractions.zeros(backend, Tv, size_eval_grid..., Nout)
    # NURBS fields not needed
    weights = nothing
    SplineGrid(
        spline_dimensions,
        control_points,
        weights,
        eval
    )
end

"""
Create a SplineGrid but with preallocated weights to define a NURBS. See The SplineGrid
constructors for more details.
"""
function NURBSGrid(spline_dimensions::NTuple{Nin, <:SplineDimension{Tv, Ti}},
        Nout::Integer
)::AbstractSplineGrid{Nin, Nout, true, Tv, Ti} where {Nin, Tv, Ti}
    nurbs_grid = SplineGrid(spline_dimensions, Nout)
    backend = get_backend(first(spline_dimensions))
    size_cp_grid = get_control_point_grid_size(spline_dimensions)
    setproperties(
        nurbs_grid; weights = KernelAbstractions.ones(backend, Tv, size_cp_grid...))
end

# Support inputting a single spline dimension instead of a tuple
function SplineGrid(spline_dimension::SplineDimension, args...; kwargs...)
    SplineGrid((spline_dimension,), args...; kwargs...)
end

function NURBSGrid(spline_dimension::SplineDimension, args...; kwargs...)
    NURBSGrid((spline_dimension,), args...; kwargs...)
end

@kernel function spline_eval_kernel(
        eval,
        @Const(basis_function_eval_all),
        @Const(sample_indices_all),
        @Const(control_points),
        @Const(weights),
        control_point_kernel_size,
        derivative_order,
        degrees
)

    # Index of the global sample point
    J = @index(Global, Cartesian)

    # Input and output dimensionality
    Nin = ndims(control_points) - 1
    Nout = size(control_points)[end]

    # Set value to 0
    for dim_out in 1:Nout
        eval[J, dim_out] = 0
    end

    is_nurbs = !isnothing(weights)
    denom = zero(eltype(eval))

    control_point_index_base = ntuple(
        dim_in -> sample_indices_all[dim_in][J[dim_in]] - degrees[dim_in] - 1, Nin)

    for I in CartesianIndices(control_point_kernel_size)
        control_point_index = ntuple(
            dim_in -> control_point_index_base[dim_in] + I[dim_in], Nin)

        # Compute basis function product
        basis_function_product = one(eltype(eval))

        for dim_in in 1:Nin
            basis_function_product *= basis_function_eval_all[dim_in][
                J[dim_in], I[dim_in], derivative_order[dim_in] + 1]
        end

        # Handle NURBS weights and denominator if is_nurbs
        basis_function_multidim = if is_nurbs
            basis_function_multidim = basis_function_product *
                                      weights[control_point_index...]
            denom += basis_function_multidim
            basis_function_multidim
        else
            basis_function_product
        end

        # Add product of multi-dimensional basis function and control point to output
        for dim_out in 1:Nout
            eval[J, dim_out] += basis_function_multidim *
                                control_points[control_point_index..., dim_out]
        end
    end

    # Divide by denominator if is_nurbs
    if is_nurbs
        for dim_out in 1:Nout
            eval[J, dim_out] /= denom
        end
    end
end

"""
    evaluate!(spline_grid::SplineGrid{Nin};
        derivative_order::NTuple{Nin, <:Integer} = ntuple(_ -> 0, Nin),
        control_points::AbstractArray = spline_grid.control_points,
        eval::AbstractArray = spline_grid.eval)

Evaluate the spline grid, that is: take the evaluated basis functions for each sample point
for each SplineDimension, and compute the output grid on each sample point combination
as a linear combination of control points with basis function products as coefficients.

If weights are supplied, compute the rational basis functions for NURBS as the control point coefficients.

Uses the `control_points` and `eval` arrays from the `spline_grid` by default,
but different arrays can be specified as a convenience for optimization algorithms.
"""
function evaluate!(spline_grid::AbstractSplineGrid{Nin, Nout, HasWeights, Tv};
        derivative_order::NTuple{Nin, <:Integer} = ntuple(_ -> 0, Nin),
        control_points::AbstractControlPointArray{Nin, Nout, Tv} = spline_grid.control_points,
        eval::AbstractArray = spline_grid.eval
)::Nothing where {Nin, Nout, HasWeights, Tv}
    (; spline_dimensions, weights) = spline_grid

    validate_partial_derivatives(spline_grid, derivative_order)
    @assert size(control_points) == size(spline_grid.control_points)
    @assert size(eval) == size(spline_grid.eval)

    basis_function_eval_all = ntuple(i -> spline_dimensions[i].eval, Nin)
    sample_indices_all = ntuple(i -> spline_dimensions[i].sample_indices, Nin)
    control_point_kernel_size = get_cp_kernel_size(spline_dimensions)
    degrees = ntuple(i -> spline_dimensions[i].degree, Nin)

    backend = get_backend(eval)
    spline_eval_kernel(backend)(
        eval,
        basis_function_eval_all,
        sample_indices_all,
        obtain(control_points),
        weights,
        control_point_kernel_size,
        derivative_order,
        degrees,
        ndrange = size(eval)[1:(end - 1)]
    )
    synchronize(backend)
    return nothing
end

@kernel function spline_eval_adjoint_kernel(
        control_points,
        @Const(basis_function_eval_all),
        @Const(sample_indices_all),
        @Const(eval),
        control_point_kernel_size,
        derivative_order,
        degrees
)

    # Index of the global sample point
    J = @index(Global, Cartesian)

    # Input and output dimensionality
    Nin = ndims(control_points) - 1
    Nout = size(control_points)[end]

    control_point_index_base = ntuple(
        dim_in -> sample_indices_all[dim_in][J[dim_in]] - degrees[dim_in] - 1, Nin)

    for I in CartesianIndices(control_point_kernel_size)
        control_point_index = ntuple(
            dim_in -> control_point_index_base[dim_in] + I[dim_in], Nin)

        # Compute basis function product
        basis_function_product = one(eltype(eval))

        for dim_in in 1:Nin
            basis_function_product *= basis_function_eval_all[dim_in][
                J[dim_in], I[dim_in], derivative_order[dim_in] + 1]
        end

        # Add product of multi-dimensional basis function and control point to output
        for dim_out in 1:Nout
            Atomix.@atomic control_points[
            control_point_index..., dim_out] += basis_function_product *
                                                eval[J, dim_out]
        end
    end
end

"""
    evaluate_adjoint!(spline_grid::AbstractSplineGrid{Nin, Nout, false, Tv};
        derivative_order::NTuple{Nin, <:Integer} = ntuple(_ -> 0, Nin),
        control_points::AbstractControlPointArray{Nin, Nout, Tv} = spline_grid.control_points,
        eval::AbstractArray = spline_grid.eval)::Nothing where {Nin, Nout, Tv}

evaluate the adjoint of the linear mapping `control_points -> eval`. This is a computation of the form
`eval -> control_points`. If we write `evaluate!(spline_grid)` as a matrix vector multiplication `eval = M * control_points`,
Then the adjoint is given by `v -> M' * v`. This mapping is used in fitting algorithms.
"""
function evaluate_adjoint!(spline_grid::AbstractSplineGrid{Nin, Nout, false, Tv};
        derivative_order::NTuple{Nin, <:Integer} = ntuple(_ -> 0, Nin),
        control_points::AbstractControlPointArray{Nin, Nout, Tv} = spline_grid.control_points,
        eval::AbstractArray = spline_grid.eval
)::Nothing where {Nin, Nout, Tv}
    @assert !is_nurbs(spline_grid) "Adjoint evaluation not supported for NURBS."
    (; spline_dimensions) = spline_grid
    validate_partial_derivatives(spline_dimensions, derivative_order)
    control_points = obtain(control_points)
    control_points .= 0
    @assert size(control_points) == size(spline_grid.control_points)
    @assert size(eval) == size(spline_grid.eval)

    basis_function_eval_all = ntuple(i -> spline_dimensions[i].eval, Nin)
    sample_indices_all = ntuple(i -> spline_dimensions[i].sample_indices, Nin)
    control_point_kernel_size = get_cp_kernel_size(spline_dimensions)
    degrees = ntuple(i -> spline_dimensions[i].degree, Nin)

    backend = get_backend(eval)
    spline_eval_adjoint_kernel(backend)(
        control_points,
        basis_function_eval_all,
        sample_indices_all,
        eval,
        control_point_kernel_size,
        derivative_order,
        degrees,
        ndrange = size(eval)[1:(end - 1)]
    )
    synchronize(backend)
    return nothing
end

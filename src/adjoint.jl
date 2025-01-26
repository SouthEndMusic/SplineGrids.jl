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

@kernel function refinement_matrix_array_mul_adjoint_kernel(
        B,
        @Const(Y),
        @Const(row_pointer_all),
        @Const(column_start_all),
        @Const(nzval_all),
        refmat_index_all
)
    # Y index
    I = @index(Global, Cartesian)

    Ndims = ndims(Y)

    column_start, n_columns = get_row_extends(
        I,
        refmat_index_all,
        row_pointer_all,
        column_start_all,
        nzval_all
    )

    for J_base in CartesianIndices(Tuple(n_columns))
        # B index
        J = ntuple(dim -> J_base[dim] + column_start[dim] - 1, Ndims)
        contrib = Y[I]
        for dim in 1:Ndims
            refmat_index = refmat_index_all[dim]
            if !iszero(refmat_index)
                row_pointer = row_pointer_all[refmat_index][I[dim]]
                contrib *= nzval_all[refmat_index][row_pointer + J_base[dim] - 1]
            end
        end
        if contrib isa Flag
            if contrib.flag
                B[J...] = contrib
            end
        else
            Atomix.@atomic B[J...] += contrib
        end
    end
end

function mult_adjoint!(
        B::AbstractArray,
        As::NTuple{N, <:RefinementMatrix},
        Y::AbstractArray,
        dims_refinement::NTuple{N, <:Integer}
) where {N}
    backend = get_backend(B)
    validate_mult_input(Y, As, B, dims_refinement)

    n_refmat = length(dims_refinement)
    refmat_index_all = ntuple(
        dim -> (dim ∈ dims_refinement) ? findfirst(==(dim), dims_refinement) : 0, ndims(Y))

    refinement_matrix_array_mul_adjoint_kernel(backend)(
        B,
        Y,
        ntuple(i -> As[i].row_pointer, n_refmat),
        ntuple(i -> As[i].column_start, n_refmat),
        ntuple(i -> As[i].nzval, n_refmat),
        refmat_index_all,
        ndrange = size(Y)
    )
    synchronize(backend)
    return nothing
end

@kernel function local_refinement_adjoint_kernel(
        refinement_values,
        @Const(control_points),
        @Const(refinement_indices)
)
    i = @index(Global, Linear)

    Nin = ndims(control_points) - 1
    Nout = size(control_points)[end]

    indices = ntuple(dim_in -> refinement_indices[i, dim_in], Nin)

    for dim_out in 1:Nout
        refinement_values[i, dim_out] = control_points[indices..., dim_out]
    end
end

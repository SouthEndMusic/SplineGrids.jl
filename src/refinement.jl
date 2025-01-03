@kernel function build_refinement_matrix_kernel(
        row_pointer,
        column_start,
        nzval,
        @Const(knots_all_old),
        knot_span_index,
        knot_new,
        degree)
    # Row index
    i = @index(Global, Linear)

    if i ≤ knot_span_index - degree
        row_pointer[i] = i
        column_start[i] = i
    elseif i ≤ knot_span_index
        α = (knot_new -
             knots_all_old[i]) /
            (knots_all_old[i + degree] -
             knots_all_old[i])
        row_pointer_ = 2i - knot_span_index + degree - 1
        row_pointer[i] = row_pointer_
        column_start[i] = i - 1
        nzval[row_pointer_] = 1 - α
        nzval[row_pointer_ + 1] = α
    else
        row_pointer_ = i + degree
        row_pointer[i] = row_pointer_
        column_start[i] = i - 1
    end
end

function RefinementMatrix(
        spline_dimension::AbstractSplineDimension{Tv},
        knot_span_index::Ti,
        knot_new
)::AbstractRefinementMatrix{Tv, Ti} where {Tv, Ti <: Integer}
    (; degree, knot_vector) = spline_dimension

    n_basis_functions = get_n_basis_functions(spline_dimension)
    n_basis_functions_new = n_basis_functions + 1
    n_nzval = n_basis_functions + degree + 1

    backend = get_backend(spline_dimension)
    row_pointer = allocate(backend, Ti, n_basis_functions_new)
    column_start = allocate(backend, Ti, n_basis_functions_new)
    nzval = KernelAbstractions.ones(backend, Tv, n_nzval)

    build_refinement_matrix_kernel(backend)(
        row_pointer,
        column_start,
        nzval,
        knot_vector.knots_all,
        knot_span_index,
        knot_new,
        degree,
        ndrange = size(row_pointer)
    )
    synchronize(backend)

    RefinementMatrix(
        n_basis_functions_new,
        n_basis_functions,
        row_pointer,
        column_start,
        nzval
    )
end

"""
    insert_knot(
    knot_vector::KnotVector, knot_new::AbstractFloat)::Tuple{KnotVector, Integer}

Create a new knot vector with the new knot of multiplicity 1.

## Inputs

  - `knot_vector`: The knot vector object to which the knot will be added
  - `knot_new`: The value of the new knot. Should not be part of the knot values already

## Outputs

  - `knot_vector_new`: The newly created knot vector with the added knot
  - `knot_span_index`: The index of the old knot vector knot span which contains the new knot
"""
function insert_knot(
        knot_vector::KnotVector, knot_new::AbstractFloat
)::Tuple{KnotVector, Integer}
    (; knot_values, multiplicities) = knot_vector
    @assert !any(knot_new .== knot_values)

    knot_values_knot_new_index = searchsortedfirst(adapt(CPU(), knot_values), knot_new)
    knot_values_new = insert(
        knot_values,
        knot_values_knot_new_index,
        knot_new
    )
    multiplicities_new = insert(
        multiplicities,
        knot_values_knot_new_index,
        1
    )

    knot_vector_new = KnotVector(
        knot_values_new,
        multiplicities_new
    )

    knot_span_index = sum(view(multiplicities, 1:(knot_values_knot_new_index - 1)))
    knot_vector_new, knot_span_index
end

"""
    insert_knot(
            spline_dimension::SplineDimension,
            knot_new::AbstractFloat;
            recompute_sample_indices::Bool = true,
            evaluate::Bool = true::Tuple{SplineDimension, RefinementMatrix}

Create a new spline dimension whose knot vector has the new knot with multiplicity 1.

## Inputs

  - `spline_dimension`: The spline dimension the new knot will be added to
  - `knot_new`: The value of the new knot
  - `recompute_sample_indices`: Whether the indices of the sample points should be recomputed after the knot insertion.
    Defaults to `true`.
  - `evaluate`: Whether the spline dimension should be evaluated after the knot insertion. Defaults to `true`.

## Outputs

  - `spline_dimension_new`: The newly created spline dimension with the same underlying memory except for the new knot vector.
  - `refinement_matrix`: The sparse matrix which expresses the basis functions from before the knot insertion in terms of
    the basis functions after the knot insertion.
"""
function insert_knot(
        spline_dimension::SplineDimension,
        knot_new::AbstractFloat;
        recompute_sample_indices::Bool = true,
        evaluate::Bool = true
)::Tuple{SplineDimension, RefinementMatrix}
    (; knot_vector) = spline_dimension

    # Add new knot to knot vector
    knot_vector_new, knot_span_index = insert_knot(knot_vector, knot_new)
    spline_dimension_new = @set spline_dimension.knot_vector = knot_vector_new

    # Compute refinement matrix
    refinement_matrix = RefinementMatrix(spline_dimension, knot_span_index, knot_new)

    # Recompute the sample indices
    if recompute_sample_indices
        set_sample_indices!(spline_dimension_new)
        evaluate && evaluate!(spline_dimension_new)
    end

    spline_dimension_new, refinement_matrix
end

"""
    insert_knot(
    spline_grid::AbstractSplineGrid{Nin, Nout},
    knot_new::AbstractFloat,
    dim_refinement::Integer;
    evaluate_spline_dimension::Bool = true,
    recompute_global_sample_indices = true
    )::Tuple{SplineGrid, RefinementMatrix} where {Nin, Nout}

Create a new spline grid where a new knot is added to the knot vector underlying the indicated spline dimension.

## Inputs

  - `spline_grid`: The spline grid to which the new knot will be added
  - `knot_new`: The value of the knot to be added
  - `dim_refinement`: The index of the spline dimension to which the knot will be added
  - `evaluate_spline_dimension`: Whether the spline dimension to which the knot is added should be evaluated.
    Defaults to `true`.
  - `recompute_global_sample_indices`: Whether the global sample indices should be recomputed after the knot insertion.
    Defaults to `true`.

## Outputs

  - `spline_grid_new`: The newly created spline grid with all the same underlying memory
    except for the updated knot vector and the control points.
  - `refinement_matrix`: The sparse matrix which expresses the basis functions from before the knot insertion in terms of
    the basis functions after the knot insertion for the knot insertion dimension.
"""
function insert_knot(
        spline_grid::AbstractSplineGrid{Nin, Nout, HasWeights, T},
        dim_refinement::Integer,
        knot_new::AbstractFloat;
        evaluate_spline_dimension::Bool = true,
        recompute_global_sample_indices = true
)::Tuple{SplineGrid, RefinementMatrix} where {Nin, Nout, HasWeights, T}
    (; spline_dimensions, control_points, sample_indices) = spline_grid

    @assert !is_nurbs(spline_grid) "Knot insertion not supported for NURBS"

    spline_dimension = spline_dimensions[dim_refinement]
    spline_dimension_new, refinement_matrix = insert_knot(
        spline_dimension,
        knot_new;
        evaluate = evaluate_spline_dimension
    )

    # Refine control points
    size_cp_grid = size(control_points)
    size_cp_grid_new = ntuple(
        dim -> (dim == dim_refinement) ? size_cp_grid[dim] + 1 :
               size_cp_grid[dim],
        ndims(control_points))

    backend = get_backend(spline_dimension)
    control_points_new = allocate(
        backend,
        T,
        size_cp_grid_new...
    )

    mul!(
        control_points_new,
        refinement_matrix,
        control_points,
        dim_refinement
    )

    spline_dimensions_new = ntuple(
        dim -> (dim == dim_refinement) ? spline_dimension_new : spline_dimensions[dim], Nin)

    if recompute_global_sample_indices
        set_global_sample_indices!(
            sample_indices,
            spline_dimensions_new,
            Nout)
    end

    spline_grid_new = @set spline_grid.control_points = control_points_new
    spline_grid_new = @set spline_grid_new.spline_dimensions = spline_dimensions_new

    spline_grid_new, refinement_matrix
end

"""
    refine(spline_dimension::SplineDimension;
    knots_new::Union{Vector{<:AbstractFloat}, Nothing} = nothing)::Tuple{SplineDimension, RefinementMatrix}

Create a new spline dimension with multiple knots added to the underlying knot vector.

## Inputs

  - `spline_dimension`: The spline dimension the new knots will be added to
  - `knots_new`: The vector of knots that will be added. Defaults to the midpoints of the knot spans of the vector underlying
    the spline dimension.

## Outputs

  - `spline_dimension_new`: The newly created spline dimension with the same underlying memory except for the new knot vector.
  - `refinement_matrix`: The sparse matrix which expresses the basis functions from before the refinement in terms of
    the basis functions after the refinement.
"""
function refine(
        spline_dimension::AbstractSplineDimension{Tv};
        knots_new::Union{Vector{<:AbstractFloat}, Nothing} = nothing
)::Tuple{SplineDimension, RefinementMatrix} where {Tv}
    (; knot_values) = spline_dimension.knot_vector
    backend = get_backend(spline_dimension)

    # Default new knots: midpoints of knot spans
    if isnothing(knots_new)
        knots_new = knot_values[1:(end - 1)] .+ diff(knot_values) / 2
    end

    # Start with refinement matrix = identity matrix
    n_basis_functions = get_n_basis_functions(spline_dimension)
    refinement_matrix = rmeye(n_basis_functions; backend, Tv)

    # Loop over knots to insert them in the knot vector and multiply the 
    # refinement matrices to obtain the refinement matrix of the whole refinement
    spline_dimension_new = spline_dimension
    for knot_new in adapt(CPU(), knots_new)
        spline_dimension_new, refinement_matrix_knot = insert_knot(
            spline_dimension_new, knot_new; recompute_sample_indices = false)
        refinement_matrix = refinement_matrix_knot * refinement_matrix
    end

    set_sample_indices!(spline_dimension_new)

    spline_dimension_new, refinement_matrix
end

"""
    refine(
    spline_grid::AbstractSplineGrid{Nin, Nout, T},
    dim_refinement::Integer;
    knots_new::Union{Vector{<:AbstractFloat}, Nothing} = nothing,
    recompute_global_sample_indices = true
    )::Tuple{SplineGrid, RefinementMatrix} where {Nin, Nout, T}

Create a new spline grid where multiple knots are added to the knot vector underlying the indicated spline dimension.

## Inputs

  - `spline_grid`: The spline grid from which one of the knot vectors will be refined
  - `dim_refinement`: The index of the spline dimension whose knot vector will be refined
  - `knots_new`: The knots that will be added. Defaults to `nothing`, which internally is translated to all midpoints of the
    non-trivial knot spans of the knot vector that will be refined.
  - `recompute_global_sample_indices`: Whether the global sample indices should be recomputed after the knot refinement.
    Defaults to `true`.

## Outputs

  - `spline_grid_new`: The newly created spline grid with all the same underlying memory
    except for the updated knot vector and the control points.
  - `refinement_matrix`: The sparse matrix which expresses the basis functions from before the knot insertion in terms of
    the basis functions after the knot insertion for the knot insertion dimension.
"""
function refine(
        spline_grid::AbstractSplineGrid{Nin, Nout, HasWeights, T},
        dim_refinement::Integer;
        knots_new::Union{Vector{<:AbstractFloat}, Nothing} = nothing,
        recompute_global_sample_indices::Bool = true
)::Tuple{SplineGrid, RefinementMatrix} where {Nin, Nout, HasWeights, T}
    (; spline_dimensions, control_points, sample_indices) = spline_grid

    @assert !is_nurbs(spline_grid) "Knot insertion not supported for NURBS"

    spline_dimension = spline_dimensions[dim_refinement]
    spline_dimension_new, refinement_matrix = refine(
        spline_dimension;
        knots_new
    )
    evaluate!(spline_dimension_new)

    # The number of control points added in the refinement dimension
    n_cp_new = size(refinement_matrix)[1] - size(refinement_matrix)[2]

    # Refine the control points
    size_cp_grid = size(control_points)
    size_cp_grid_new = ntuple(
        dim -> (dim == dim_refinement) ? size_cp_grid[dim] + n_cp_new :
               size_cp_grid[dim],
        ndims(control_points))

    backend = get_backend(spline_dimension)
    control_points_new = allocate(
        backend,
        T,
        size_cp_grid_new...
    )

    mul!(
        control_points_new,
        refinement_matrix,
        control_points,
        dim_refinement
    )

    spline_dimensions_new = ntuple(
        dim -> (dim == dim_refinement) ? spline_dimension_new : spline_dimensions[dim], Nin)

    if recompute_global_sample_indices
        set_global_sample_indices!(
            sample_indices,
            spline_dimensions_new,
            Nout
        )
    end

    spline_grid_new = @set spline_grid.control_points = control_points_new
    spline_grid_new = @set spline_grid_new.spline_dimensions = spline_dimensions_new

    spline_grid_new, refinement_matrix
end
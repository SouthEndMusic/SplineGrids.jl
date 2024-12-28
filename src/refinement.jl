# The new knot lies between old knots (knot_new_index - 1) and knot_new index
# The old basis functions that have support there are (knot_new_index - degree - 1), ..., (knot_new_index - 1)
function build_refinement_matrix(
        spline_dimension::SplineDimension, knot_new_index::Integer)::SparseMatrixCSC
    (; degree, knot_vector) = spline_dimension
    (; knots_all) = knot_vector

    knot_new = knots_all[knot_new_index]
    n_basis_functions = get_n_basis_functions(spline_dimension)
    refinement_matrix = spzeros(n_basis_functions, n_basis_functions - 1)

    for i in 1:(n_basis_functions - 1)
        if i ≤ knot_new_index - degree - 2
            # Basis functions left of the knot insertion that remain the same
            refinement_matrix[i, i] = 1
        elseif i ≥ knot_new_index
            # Basis functions right of the knot insertion that remain the same
            refinement_matrix[i + 1, i] = 1
        else
            # Basis functions whose support contains the new knot
            refinement_matrix[i, i] = (knot_new -
                                       knots_all[i]) /
                                      (knots_all[i + degree + 1] -
                                       knots_all[i])
            refinement_matrix[i + 1, i] = (knots_all[i + degree + 2] -
                                           knot_new) /
                                          (knots_all[i + degree + 2] -
                                           knots_all[i + 1])
        end
    end

    refinement_matrix
end

function refine_control_points!(
        control_points_new::AbstractArray,
        control_points::AbstractArray,
        refinement_matrix::AbstractMatrix,
        dim_refinement::Integer
)::Nothing
    size_control_points = size(control_points)

    # The refinement matrix has to be applied for every index in every dimension apart from
    # the refinement dimension
    refinement_indices = ntuple(
        dim_in -> dim_in == dim_refinement ? 1 : size_control_points[dim_in], ndims(control_points))

    for I in CartesianIndices(refinement_indices)
        indices = ntuple(
            dim_in -> dim_in == dim_refinement ? Colon() : I[dim_in], ndims(control_points))
        mul!(
            view(control_points_new, indices...),
            refinement_matrix,
            view(control_points, indices...)
        )
    end
end

"""
    insert_knot!(
    knot_vector::KnotVector, knot_new::AbstractFloat)::Integer

Add a new knot to the knot vector with multiplicity 1.

## Inputs

  - `knot_vector`: The knot vector object to which the knot will be added
  - `knot_new`: The value of the new knot. Should not be part of the knot values already

## Outputs

  - `knot_new_index`: The index of the new knot in `knot_vector.knots_all`
"""
function insert_knot!(
        knot_vector::KnotVector, knot_new::AbstractFloat)::Integer
    (; knot_values, multiplicities, knots_all) = knot_vector
    @assert knot_new ∉ knot_values

    # Index of the new knot in the knot values (without multiplicities)
    knot_new_index_knot_values = searchsortedfirst(knot_values, knot_new)
    insert!(
        knot_values,
        knot_new_index_knot_values,
        knot_new
    )
    insert!(
        multiplicities,
        knot_new_index_knot_values,
        1
    )
    # Index of the new knot in all knots (with multiplicities)
    knot_new_index = searchsortedfirst(knots_all, knot_new)
    insert!(
        knots_all,
        knot_new_index,
        knot_new
    )
    return knot_new_index
end

"""
    insert_knot!(
            spline_dimension::SplineDimension,
            knot_new::AbstractFloat;
            recompute_sample_indices::Bool = true,
            evaluate::Bool = true)::SparseMatrixCSC

Add a new knot to the knot vector underlying the spline dimension.

## Inputs

  - `spline_dimension`: The spline dimension the new knot will be added to
  - `knot_new`: The value of the new knot
  - `recompute_sample_indices`: Whether the indices of the sample points should be recomputed after the knot insertion.
    Defaults to `true`.
  - `evaluate`: Whether the spline dimension should be evaluated after the knot insertion. Defaults to `true`.

## Outputs

  - `refinement_matrix`: The sparse matrix which expresses the basis functions from before the knot insertion in terms of
    the basis functions after the knot insertion.
"""
function insert_knot!(
        spline_dimension::SplineDimension,
        knot_new::AbstractFloat;
        recompute_sample_indices::Bool = true,
        evaluate::Bool = true
)::SparseMatrixCSC
    (; knot_vector) = spline_dimension

    # Add new knot to knot vector
    knot_new_index = insert_knot!(knot_vector, knot_new)

    # Compute refinement matrix
    refinement_matrix = build_refinement_matrix(spline_dimension, knot_new_index)

    # Recompute the sample indices
    if recompute_sample_indices
        reset_sample_indices!(spline_dimension)
        evaluate && evaluate!(spline_dimension)
    end

    refinement_matrix
end

"""
    insert_knot!(
    spline_grid::AbstractSplineGrid{Nin, Nout},
    knot_new::AbstractFloat,
    dim_refinement::Integer;
    evaluate_spline_dimension::Bool = true,
    recompute_global_sample_indices = true
    )::Tuple{AbstractSplineGrid, SparseMatrixCSC} where {Nin, Nout}

Add a new knot to the knot vector underlying the indicated spline dimension.
NOTE: A new `SplineGrid` object is created, with the same memory for all underlying data
except for the control points.

## Inputs

  - `spline_grid`: The spline grid to which the new knot will be added
  - `knot_new`: The value of the knot to be added
  - `dim_refinement`: The index of the spline dimension to which the knot will be added
  - `evaluate_spline_dimension`: Whether the spline dimension to which the knot is added should be evaluated.
    Defaults to `true`.
  - `recompute_global_sample_indices`: Whether the global sample indices should be recomputed after the knot insertion.
    Defaults to `true`.

## Outputs

  - `spline_grid_new`: The newly created spline grid with all the same underlying memory except for the control points.
  - `refinement_matrix`: The sparse matrix which expresses the basis functions from before the knot insertion in terms of
    the basis functions after the knot insertion for the knot insertion dimension.
"""
function insert_knot!(
        spline_grid::AbstractSplineGrid{Nin, Nout},
        dim_refinement::Integer,
        knot_new::AbstractFloat;
        evaluate_spline_dimension::Bool = true,
        recompute_global_sample_indices = true
)::Tuple{AbstractSplineGrid, SparseMatrixCSC} where {Nin, Nout}
    (; spline_dimensions, control_points, sample_indices) = spline_grid

    @assert !is_nurbs(spline_grid) "Knot insertion not supported for NURBS"

    spline_dimension = spline_dimensions[dim_refinement]
    refinement_matrix = insert_knot!(
        spline_dimension,
        knot_new;
        evaluate = evaluate_spline_dimension)

    # Refine control points
    size_cp_grid = size(control_points)
    size_cp_grid_new = ntuple(
        dim_in -> dim_in == dim_refinement ? size_cp_grid[dim_in] + 1 :
                  size_cp_grid[dim_in],
        ndims(control_points))

    control_points_new = zeros(size_cp_grid_new)

    refine_control_points!(
        control_points_new,
        control_points,
        refinement_matrix,
        dim_refinement)

    if recompute_global_sample_indices
        set_global_sample_indices!(
            sample_indices,
            spline_dimensions,
            control_points_new)
    end

    spline_grid_new = SplineGrid(
        spline_dimensions,
        control_points_new,
        spline_grid.eval,
        spline_grid.sample_indices,
        spline_grid.basis_function_products
    )
    spline_grid_new, refinement_matrix
end

"""
    refine!(spline_dimension::SplineDimension;
    knots_new::Union{Vector{<:AbstractFloat}, Nothing} = nothing)::SparseMatrixCSC

Add multiple knots at once to the knot vector underlying the spline dimension.

## Inputs

  - `spline_dimension`: The spline dimension the new knots will be added to
  - `knots_new`: The vector of knots that will be added. Defaults to the midpoints of the knot spans of the vector underlying
    the spline dimension.

## Outputs

  - `refinement_matrix`: The sparse matrix which expresses the basis functions from before the refinement in terms of
    the basis functions after the refinement.
"""
function refine!(
        spline_dimension::SplineDimension;
        knots_new::Union{Vector{<:AbstractFloat}, Nothing} = nothing
)::SparseMatrixCSC
    (; knot_values) = spline_dimension.knot_vector

    # Default new knots: midpoints of knot spans
    if isnothing(knots_new)
        knots_new = knot_values[1:(end - 1)] .+ diff(knot_values) / 2
    end

    refinement_matrix = one(eltype(knots_new))

    # Loop over knots to insert them in the knot vector and multiply the 
    # refinement matrices to obtain the refinement matrix of the whole refinement
    for knot_new in knots_new
        refinement_matrix_knot = insert_knot!(
            spline_dimension, knot_new; recompute_sample_indices = false)
        refinement_matrix = refinement_matrix_knot * refinement_matrix
    end

    reset_sample_indices!(spline_dimension)

    refinement_matrix
end

"""
    refine!(
    spline_grid::AbstractSplineGrid{Nin, Nout},
    dim_refinement::Integer;
    knots_new::Union{Vector{<:AbstractFloat}, Nothing} = nothing,
    recompute_global_sample_indices = true
    )::Tuple{AbstractSplineGrid, SparseMatrixCSC} where {Nin, Nout}

Add multiple knots to the knot vector underlying the indicated spline dimension.
NOTE: A new `SplineGrid` object is created, with the same memory for all underlying data
except for the control points.

## Inputs

  - `spline_grid`: The spline grid from which one of the knot vectors will be refined
  - `dim_refinement`: The index of the spline dimension whose knot vector will be refined
  - `knots_new`: The knots that will be added. Defaults to `nothing`, which internally is translated to all midpoints of the
    non-trivial knot spans of the knot vector that will be refined.
  - `recompute_global_sample_indices`: Whether the global sample indices should be recomputed after the knot refinement.
    Defaults to `true`.

## Outputs

  - `spline_grid_new`: The newly created spline grid with all the same underlying memory except for the control points.
  - `refinement_matrix`: The sparse matrix which expresses the basis functions from before the knot insertion in terms of
    the basis functions after the knot insertion for the knot insertion dimension.
"""
function refine!(
        spline_grid::AbstractSplineGrid,
        dim_refinement::Integer;
        knots_new::Union{Vector{<:AbstractFloat}, Nothing} = nothing,
        recompute_global_sample_indices::Bool = true
)::Tuple{AbstractSplineGrid, SparseMatrixCSC}
    (; spline_dimensions, control_points, sample_indices) = spline_grid
    spline_dimension = spline_dimensions[dim_refinement]
    refinement_matrix = refine!(spline_dimension; knots_new)
    evaluate!(spline_dimension)

    # The number of control points added in the refinement dimension
    n_cp_new = size(refinement_matrix)[1] - size(refinement_matrix)[2]

    # Refine the control points
    size_cp_grid = size(control_points)
    size_cp_grid_new = ntuple(
        dim_in -> dim_in == dim_refinement ? size_cp_grid[dim_in] + n_cp_new :
                  size_cp_grid[dim_in],
        ndims(control_points))
    control_points_new = zeros(size_cp_grid_new)
    refine_control_points!(
        control_points_new,
        control_points,
        refinement_matrix,
        dim_refinement
    )

    spline_grid_new = SplineGrid(
        spline_dimensions,
        control_points_new,
        spline_grid.eval,
        spline_grid.sample_indices,
        spline_grid.basis_function_products
    )

    if recompute_global_sample_indices
        set_global_sample_indices!(sample_indices, spline_dimensions, control_points_new)
    end

    spline_grid_new, refinement_matrix
end
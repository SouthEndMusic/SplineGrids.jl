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
            refinement_matrix[i, i] = 1
        elseif i ≥ knot_new_index
            refinement_matrix[i + 1, i] = 1
        else
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

"""
    insert_knot!(
            spline_dimension::SplineDimension,
            knot_new::AbstractFloat;
            recompute_sample_indices = true)::SparseMatrixCSC

Add a new knot to the knot vector underlying the spline dimension.

## Inputs

  - `spline_dimension`: The spline dimension the new knot will be added to
  - `knot_new`: The value of the new knot
  - `recompute_sample_indices`: Whether the indices of the sample points should be recomputed after the knot insertion.
    Defaults to `true`.

## Outputs

  - `refinement_matrix`: The sparse matrix which expresses the basis functions from before the knot insertion in terms of
    the basis functions after the knot insertion.
"""
function insert_knot!(
        spline_dimension::SplineDimension,
        knot_new::AbstractFloat;
        recompute_sample_indices = true
)::SparseMatrixCSC
    (; knot_vector) = spline_dimension

    # Add new knot to knot vector
    knot_new_index = insert_knot!(knot_vector, knot_new)

    # Compute refinement matrix
    refinement_matrix = build_refinement_matrix(spline_dimension, knot_new_index)

    # Recompute the sample indices
    if recompute_sample_indices
        reset_sample_indices!(spline_dimension)
    end

    refinement_matrix
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
function refine!(spline_dimension::SplineDimension;
        knots_new::Union{Vector{<:AbstractFloat}, Nothing} = nothing)::SparseMatrixCSC
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
"""
    SplineDimension(degree, knot_vector, sample_points, sample_indices)

Defines the set of basis functions for a single dimension, and how it is sampled.

## Arguments

  `degree`: The degree of the piecewise polynomial basis functions.
  `knot_vector`: The knot vector on which the basis functions are defined.
  `sample_points`: The points in the domain of the basis functions where they are sampled. Must
     lie within the boundaries of the knot vector.
  `sample_indices`: The indices `i` of the sample points `t`` in the knot vector such that `knot_vector.knots[i] ≤ t < knot_vector.knots[i + 1]``
  `eval`: A matrix of shape `(length(sample_points), degree + 1)`, with per sample point the values of those basis functions 
     whose support the sample point is in. 
"""
struct SplineDimension{K, M, S <: AbstractVector, I <: AbstractVector{<:Integer}, E <: AbstractMatrix}
    degree::Int
    knot_vector::KnotVector{K, M}
    sample_points::S
    sample_indices::I
    eval::E
end

get_index(knot_vector::KnotVector, t, d) = clamp(searchsortedlast(knot_vector.knots_all, t), 1, length(knot_vector.knots_all) - d - 1)

function SplineDimension(n::Integer, d::Integer, N::Integer)::SplineDimension
    knot_vector = KnotVector(n, d)
    (; knot_values) = knot_vector
    sample_points = range(first(knot_values), last(knot_values); length = N)
    sample_indices = get_index.(Ref(knot_vector), sample_points, d)
    eval = zeros(N, d + 1)
    s = SplineDimension(d, knot_vector, sample_points, sample_indices, eval)
    evaluate!(s)
    s
end

"""
    evaluate!(spline_dimension)

Per sample point, get the value of the `spline_dimension.degree + 1` basis functions that have a
non-zero value for that sample point. This is based on the Cox-de Boor recursion formula.

The l-th sample point `t` has sample index `i`, meaning that `t ∈ [tᵢ, tᵢ₊₁)`.
Therefore `Bᵢ₀(t) = 1, Bⱼ₀(t) = 0 for j ≠ i`.
For degree `k`, `t` is in the domain of `Bⱼₖ` which is `[tⱼ, tⱼ₊ₖ₊₁)`, for `j = i - k, ..., i`.

## Arguments

- `spline_dimension`
"""
function evaluate!(spline_dimension::SplineDimension)::Nothing
    (; degree, knot_vector, sample_points, sample_indices, eval) = spline_dimension
    (; knots_all) = knot_vector

    # TODO: Preallocate or get rid of entirely
    B_prev = zeros(degree + 1)

    for (l, (t, i)) in enumerate(zip(sample_points, sample_indices))
        eval[l, 1] = 1
        for k in 1:degree
            B_prev .= eval[l, :]
            eval[l, :] .= 0
            for k_ in 1:k
                B_old = B_prev[k_]
                t_min = knots_all[i + k_ - k]
                t_max = knots_all[i + k_]
                Δt = t_max - t_min
                # Additions sum to B_old => partition of unity
                frac = B_old / Δt
                eval[l, k_] += frac * (t_max - t)
                eval[l, k_ + 1] += frac * (t - t_min)
            end
        end
    end
    return nothing
end

function decompress(spline_dimension)
    (; sample_indices, degree, knot_vector, eval) = spline_dimension
    n_sample_points = length(sample_indices)
    n_knots = sum(knot_vector.multiplicities)
    n_basis_functions = n_knots - degree - 1
    out = zeros(n_sample_points, n_basis_functions)

    for (l,i) in enumerate(sample_indices)
        out[l, (i-degree):i] .= eval[l, :]
    end

    out
end
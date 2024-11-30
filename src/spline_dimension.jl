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
struct SplineDimension{
    K, M, S <: AbstractVector, I <: AbstractVector{<:Integer}, E <: AbstractMatrix}
    degree::Int
    knot_vector::KnotVector{K, M}
    sample_points::S
    sample_indices::I
    eval::E
end

function get_index(knot_vector::KnotVector, t, d)
    clamp(searchsortedlast(knot_vector.knots_all, t),
        1, length(knot_vector.knots_all) - d - 1)
end

"""
    SplineDimension(n_basis_functions::Integer, degree::Integer, n_sample_points::Integer; kwargs...)::SplineDimension

Constructor for a SplineDimension. For now the sample points are equispaced on the extent of the knot vector.
Key word arguments are passed to the KnotVector constructor.
"""
function SplineDimension(n_basis_functions::Integer, degree::Integer, n_sample_points::Integer; kwargs...)::SplineDimension
    knot_vector = KnotVector(n_basis_functions, degree; kwargs...)
    (; knot_values) = knot_vector
    sample_points = range(first(knot_values), last(knot_values); length = n_sample_points)
    sample_indices = get_index.(Ref(knot_vector), sample_points, degree)
    eval = zeros(n_sample_points, degree + 1)
    s = SplineDimension(degree, knot_vector, sample_points, sample_indices, eval)
    evaluate!(s)
    s
end

@kernel function spline_dimension_kernel!(
        eval, @Const(knots_all), @Const(sample_points), @Const(sample_indices), degree)
    l = @index(Global, Linear)
    t = sample_points[l]
    i = sample_indices[l]

    eval[l, 1] = one(eltype(eval))
    for k in 1:degree
        @print()
        B_old = eval[l, 1]
        eval[l, 1] = zero(eltype(eval))
        for k_ in 1:k
            t_min = knots_all[i + k_ - k]
            t_max = knots_all[i + k_]
            Δt = t_max - t_min
            frac = B_old / Δt
            B_old = eval[l, k_ + 1] # Value for next iteration
            # Additions sum to B_old => partition of unity
            eval[l, k_] += frac * (t_max - t)
            eval[l, k_ + 1] = frac * (t - t_min)
        end
    end
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
    n_samples = (length(sample_points),)

    backend = get_backend(eval)
    spline_dimension_kernel!(backend)(
        eval, knots_all, sample_points, sample_indices,
        degree, ndrange = n_samples)
    synchronize(backend)
    return nothing
end

"""
Transform `spline_dimension.eval` into a matrix of shape `(n_sample_points, n_points - degree - 1)`
which explicitly gives the value for each basis function at each sample point.
"""
function decompress(spline_dimension)
    (; sample_indices, degree, knot_vector, eval) = spline_dimension
    n_sample_points = length(sample_indices)
    n_knots = sum(knot_vector.multiplicities)
    n_basis_functions = n_knots - degree - 1
    out = zeros(n_sample_points, n_basis_functions)

    for (l, i) in enumerate(sample_indices)
        out[l, (i - degree):i] .= eval[l, :]
    end

    out
end
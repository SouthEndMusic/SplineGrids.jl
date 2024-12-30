"""
    SplineDimension(
        degree, 
        max_derivative_order, 
        knot_vector, 
        sample_points, 
        sample_indices, 
        eval, 
        eval_prev)

Defines the set of basis functions for a single dimension, and how it is sampled.

## Arguments

  - `degree`: The degree of the piecewise polynomial basis functions.
  - `max_derivative_order`: The maximum derivative order of the basis functions that will be computed.
  - `knot_vector`: The knot vector on which the basis functions are defined.
  - `sample_points`: The points in the domain of the basis functions where they are sampled. Must
  - lie within the boundaries of the knot vector.
  - `sample_indices`: The indices `i` of the sample points `t` in the knot vector such that `knot_vector.knots[i] ≤ t < knot_vector.knots[i + 1]``
  - `eval`: An array of shape `(length(sample_points), degree + 1, max_derivative + 1)`, with per sample point the values of those basis functions
  - `eval_prev`: Helper array for intermediate results in the basis function computations
    whose support the sample point is in, and the derivatives if requested.
"""
struct SplineDimension{
    K <: AbstractKnotVector{T} where {T},
    S <: AbstractVector{T} where {T},
    I <: AbstractVector{<:Integer},
    E <: AbstractArray{T, 3} where {T},
    T <: AbstractFloat
} <: AbstractSplineDimension{T}
    degree::Int
    max_derivative_order::Int
    knot_vector::K
    sample_points::S
    sample_indices::I
    eval::E
    eval_prev::E
    function SplineDimension(
            degree,
            max_derivative_order,
            knot_vector,
            sample_points,
            sample_indices,
            eval,
            eval_prev
    )
        new{
            typeof(knot_vector),
            typeof(sample_points),
            typeof(sample_indices),
            typeof(eval),
            eltype(eval)
        }(
            degree,
            max_derivative_order,
            knot_vector,
            sample_points,
            sample_indices,
            eval,
            eval_prev
        )
    end
end

"""
    SplineDimension(
        n_basis_functions::Integer,
        degree::Integer,
        n_sample_points::Integer;
        max_derivative_order::Integer = 0,
        knot_vector::Union{Nothing, KnotVector} = nothing,
        backend::Backend = CPU(),
        float_type::Type{T} = Float32,
        kwargs...)

Constructor for a SplineDimension. Optionally a `knot_vector` kwarg can be passed, otherwise a default knot vector is generated.
For now by default the sample points are evenly spaced on the extent of the knot vector.
Key word arguments are passed to the KnotVector constructor.

## Inputs

  - `n_basis_functions`: The number of basis functions of this spline dimension
  - `degree`: The degree of the basis functions of this spline dimension
  - `n_sample_points`: The number of points at which the domain of the basis functions will be sampled
  - `max_derivative_order`: The maximum derivative order of the basis functions that will be computed in `evaluate!`.
    Defaults to `0`.
    `knot_vector`: A knot vector on which the basis functions will be defined. Defaults to `nothing`, which means that
    a default clamped/open equally spaced knot vector will be defined.
  - `backend`: The KernelAbstractions backend of the arrays in the object. Defaults to `CPU()`.
    NOTE: If a knot vector is supplied, its backend takes precedence.
  - `float_type`: The type of all floating point arrays. Defaults to `Float32`.
"""
function SplineDimension(
        n_basis_functions::Integer,
        degree::Integer,
        n_sample_points::Integer;
        max_derivative_order::Integer = 0,
        knot_vector::Union{Nothing, KnotVector} = nothing,
        backend::Backend = CPU(),
        float_type::Type{T} = Float32,
        kwargs...
)::SplineDimension where {T <: AbstractFloat}
    @assert 0≤max_derivative_order≤degree "The max_degree must be positive and derivatives order higher than `degree` are all 0."
    if isnothing(knot_vector)
        knot_vector = KnotVector(
            n_basis_functions,
            degree;
            backend,
            float_type,
            kwargs...)
    else
        @assert length(knot_vector.knots_all)==n_basis_functions + degree + 1 "Incompatible knot vector supplied."
        backend = get_backend(knot_vector.knot_values)
    end

    sample_points = adapt(
        backend,
        float_type.(
            range(knot_vector.extent...;
            length = n_sample_points
        )
        )
    )
    sample_indices = allocate(
        backend,
        Int,
        n_sample_points
    )

    eval = allocate(
        backend,
        float_type,
        n_sample_points,
        degree + 1,
        max_derivative_order + 1
    )
    eval_prev = similar(eval)

    s = SplineDimension(
        degree,
        max_derivative_order,
        knot_vector,
        sample_points,
        sample_indices,
        eval,
        eval_prev
    )
    set_sample_indices!(s)
    evaluate!(s)
    s
end

# NOTE: This kernel can be further optimized; the latter 2
# loops over CartesianIndices can be restricted
@kernel function spline_dimension_kernel(
        eval,
        eval_prev,
        @Const(knots_all),
        @Const(sample_points),
        @Const(sample_indices),
        degree,
        max_derivative_order
)
    l = @index(Global, Linear)
    t = sample_points[l]
    i = sample_indices[l]

    # Clear the eval_prev array for this sample point
    for I in CartesianIndices(size(eval_prev)[2:3])
        eval_prev[l, I] = zero(eltype(eval))
    end

    # Degree 0 basis function value
    eval[l, 1, 1] = one(eltype(eval_prev))
    eval_prev[l, 1, 1] = one(eltype(eval_prev))

    # Loop over successive degrees
    for k in 1:degree
        # Clear the eval array for this sample point
        for I in CartesianIndices(size(eval)[2:3])
            eval[l, I] = zero(eltype(eval))
        end
        # Loop over successive basis functions
        for k_ in 1:k
            t_min = knots_all[i + k_ - k]
            t_max = knots_all[i + k_]
            Δt = t_max - t_min
            B_prev = eval_prev[l, k_, 1]
            frac = B_prev / Δt
            # Additions sum to B_prev => partition of unity
            eval[l, k_, 1] += frac * (t_max - t)
            eval[l, k_ + 1, 1] = frac * (t - t_min)

            # Compute derivatives
            for derivative_order in 1:(max_derivative_order + k - degree)
                B_prev = eval_prev[l, k_, derivative_order]
                deriv_contribution = B_prev * k / Δt
                eval[l, k_, derivative_order + 1] -= deriv_contribution
                eval[l, k_ + 1, derivative_order + 1] = deriv_contribution
            end
        end

        # Copy from eval to eval_prev if this isn't the last k
        if k != degree
            for I in CartesianIndices(size(eval_prev)[2:3])
                eval_prev[l, I] = eval[l, I]
            end
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
    (; degree, max_derivative_order, knot_vector, sample_points, sample_indices, eval, eval_prev) = spline_dimension
    (; knots_all) = knot_vector
    n_samples = (length(sample_points),)

    backend = get_backend(eval)
    spline_dimension_kernel(backend)(
        eval, eval_prev, knots_all, sample_points, sample_indices,
        degree, ndrange = n_samples, max_derivative_order)
    synchronize(backend)
    return nothing
end

"""
Transform `spline_dimension.eval` into a matrix of shape `(n_sample_points, n_points - degree - 1)`
which explicitly gives the value for each basis function at each sample point.
"""
function decompress(
        spline_dimension::AbstractSplineDimension{T}; derivative_order = 0) where {T}
    backend = get_backend(spline_dimension)

    (; sample_indices, degree, eval, max_derivative_order) = spline_dimension
    @assert derivative_order ≤ max_derivative_order
    n_sample_points = length(sample_indices)
    n_basis_functions = get_n_basis_functions(spline_dimension)
    out = KernelAbstractions.zeros(backend, T, n_sample_points, n_basis_functions)

    for (l, i) in enumerate(sample_indices)
        out[l, (i - degree):i] .= eval[l, :, derivative_order + 1]
    end

    out
end

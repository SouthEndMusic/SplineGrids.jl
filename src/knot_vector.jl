"""
    KnotVector(knots, multiplicities)

Defines a knot vector.

## Arguments

  - `knot_values`: The values in the knot vector. Must be strictly increasing.
  - `multiplicities`: The multiplicity of each knot in `knots`.
  - `knots_all`: The explicit knot vector with all knot repeats.
"""
struct KnotVector{K <: AbstractVector{<:Number}, M <: AbstractVector{<: Integer}}
    knot_values::K
    multiplicities::M
    knots_all::K
end

function KnotVector(n::Integer, d::Integer)::KnotVector
    # For now create a clamped knot vector on [0, 1].
    knot_values = cumsum(rand(n))
    knot_values .-= knot_values[1]
    knot_values ./= knot_values[end] - knot_values[1]
    knots_all = zeros(d)
    append!(knots_all, knot_values)
    append!(knots_all, ones(d))
    multiplicities = ones(Int, n)
    multiplicities[1] = d + 1
    multiplicities[end] = d + 1
    KnotVector(knot_values, multiplicities, knots_all)
end
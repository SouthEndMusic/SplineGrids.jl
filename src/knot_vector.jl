"""
    KnotVector(knots, multiplicities)

Defines a knot vector.

## Arguments

  - `knot_values`: The values in the knot vector. Must be strictly increasing.
  - `multiplicities`: The multiplicity of each knot in `knots`.
"""
struct KnotVector{K <: AbstractVector{<:Number}, M <: AbstractVector{<:Integer}}
    knot_values::K
    multiplicities::M
    knots_all::K
    function KnotVector(knot_values, multiplicities)
        @assert length(knot_values)==length(multiplicities) "knot_values and multiplicities must be of the same length."
        @assert knot_values==sort(knot_values) "knot_values must be sorted."
        @assert allunique(knot_values) "knot_values must be unique."

        # Construct explicit knot vector
        knots_all = similar(knot_values, sum(multiplicities))
        i = 1
        for (knot_value, multiplicity) in zip(knot_values, multiplicities)
            knots_all[i:(i + multiplicity - 1)] .= knot_value
            i += multiplicity
        end
        new{typeof(knot_values), typeof(multiplicities)}(
            knot_values, multiplicities, knots_all)
    end
end

"""
    KnotVector(n_basis_functions::Integer, degree::Integer; extent::Tuple{Number, Number} = (0,1), distribution::Symbol = :equispaced)

Construct a clamped knot vector, i.e. the multiplicity of the first and last knot is degree + 1 and the
other multiplicities are 1.

## Arguments

  - `n_basis_functions`: The number of basis functions that will be defined on this knot vector
  - `degree`: The degree of the basis functions that will be defined on this knot vector

## Keyword Arguments

  - `extent`: A tuple (t_min, t_max) defining the extend of the knot vector
  - `distribution`: The distribution of the internal knots. The options are :equispaced or :random
"""
function KnotVector(
        n_basis_functions::Integer, degree::Integer; extent::Tuple{Number, Number} = (0, 1),
        distribution::Symbol = :equispaced)::KnotVector
    @assert n_basis_functions - degree ≥ 1
    n_knot_values = n_basis_functions - degree + 1
    if distribution == :random
        knot_values = cumsum(rand(n_knot_values))
        knot_values .-= knot_values[1]
        knot_values ./= knot_values[end] - knot_values[1]
        knot_values .*= extent[2] - extent[1]
        knot_values .+= extent[1]
    elseif distribution == :equispaced
        knot_values = collect(range(extent..., length = n_knot_values))
    else
        error("Unsupported knot distribution type $distribution.")
    end

    # For now always create a clamped knot vector.
    multiplicities = ones(Int, n_knot_values)
    multiplicities[1] = degree + 1
    multiplicities[end] = degree + 1
    KnotVector(knot_values, multiplicities)
end

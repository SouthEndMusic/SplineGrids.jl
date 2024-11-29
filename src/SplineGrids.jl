module SplineGrids

"""
    KnotVector(knots, multiplicities)

Defines a knot vector.

## Arguments

  - `knots`: The values in the knot vector. Must be strictly increasing.
  - `multiplicities`: The multiplicity of each knot in `knots`.
"""
struct KnotVector{K <: AbstractVector{<:Number}, M <: AbstractVector{<:Integer}}
    knots::K
    multiplicities::M
end

function KnotVector(n::Integer, d::Integer)::KnotVector
    # For now create a clamped knot vector
    knots = range(1.0, n, length = n)
    multiplicities = ones(Int, n)
    multiplicities[1] = d
    multiplicities[end] = d
    KnotVector(knots, multiplicities)
end

"""
    SplineDimension(degree, knot_vector, sample_points, sample_indices)

Defines the set of basis functions for a single dimension, and how it is sampled.

## Arguments

  `degree`: The degree of the piecewise polynomial basis functions.
  `knot_vector`: The knot vector on which the basis functions are defined.
  `sample_points`: The points in the domain of the basis functions where they are sampled. Must
     lie within the boundaries of the knot vector.
  `sample_indices`: The indices i of the sample points t in the knot vector such that knot_vector[i] ≤ t < knot_vector[i + 1]
  `eval`: A matrix of shape (length(sample_points), degree + 1), with per sample point the values of those basis functions 
     whose support the sample point is in. 
"""
struct SplineDimension{K, M, S <: AbstractVector, I <: AbstractVector, E <: AbstractMatrix}
    degree::Int
    knot_vector::KnotVector{K, M}
    sample_points::S
    sample_indices::I
    eval::E
end

get_index(v::AbstractVector, x) = clamp(searchsortedfirst(v, x) - 1, 1, length(v) - 1)

function SplineDimension(n::Integer, d::Integer, N::Integer)::SplineDimension
    knot_vector = KnotVector(n, d)
    (; knots) = knot_vector
    sample_points = range(first(knots), last(knots); length = N)
    sample_indices = get_index.(Ref(knots), sample_points)
    eval = zeros(n, d + 1)
    s = SplineDimension(d, knot_vector, sample_points, sample_indices, eval)
    evaluate!(s)
    s
end

function evaluate!(spline_dimension::SplineDimension)::Nothing
    (; degree, knot_vector, sample_points, sample_indices, eval) = spline_dimension

    for (i, t) in enumerate(sample_points)
        eval[i, 1] = 1
        for d in 1:(degree + 1)

        end
    end
    return nothing
end

struct SplineGrid{S <: SplineDimension, C <: AbstractArray, W <: Union{AbstractArray, Nothing}, E <: AbstractArray}
    spline_dimensions::Vector{S}
    control_points::C
    weights::W
    eval::E
    function SplineGrid(spline_dimensions, control_points, weights, eval)
        # TODO: Add validation of combination of control points, weights, and basis functions
        new{eltype(spline_dimensions), typeof(control_points), typeof(weights), typeof(eval)}
    end
end

# TODO: Add indexing functionality such that indexing a spline_grid implies indexing spline_grid.eval

function SplineGrid(spline_dimensions::Vector{<:SplineDimension}, control_points::AbstractVector)::SplineGrid
    # TODO: allocate eval array
    weights = nothing
    SplineGrid(spline_dimensions, control_points, weights, eval)
end

function evaluate!(spline_grid::SplineGrid)::Nothing
    # TODO
    return nothing
end

export KnotVector, SplineDimension, SplineGrid

end # module SplineGrids
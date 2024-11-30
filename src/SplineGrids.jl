module SplineGrids

using KernelAbstractions
using RecipesBase

include("knot_vector.jl")
include("spline_dimension.jl")
include("plot_rec.jl")

struct SplineGrid{S <: SplineDimension, C <: AbstractArray,
    W <: Union{AbstractArray, Nothing}, E <: AbstractArray}
    spline_dimensions::Vector{S}
    control_points::C
    weights::W
    eval::E
    function SplineGrid(spline_dimensions, control_points, weights, eval)
        # TODO: Add validation of combination of control points, weights, and basis functions
        new{eltype(spline_dimensions), typeof(control_points),
            typeof(weights), typeof(eval)}
    end
end

# TODO: Add indexing functionality such that indexing a spline_grid implies indexing spline_grid.eval

function SplineGrid(spline_dimensions::Vector{<:SplineDimension},
        control_points::AbstractVector)::SplineGrid
    # TODO: allocate eval array
    weights = nothing
    SplineGrid(spline_dimensions, control_points, weights, eval)
end

function evaluate!(spline_grid::SplineGrid)::Nothing
    # TODO
    return nothing
end

export KnotVector, SplineDimension, SplineGrid, evaluate!

end # module SplineGrids

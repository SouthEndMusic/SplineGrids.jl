struct SplineGrid{S <: SplineDimension, C <: AbstractArray, W <: Union{AbstractArray, Nothing}, E <: AbstractArray}
    spline_dimensions::Vector{S}
    control_points::C
    weights::W
    eval::E
    function SplineGrid(spline_dimensions, control_points, weights, eval)
        # TODO: Add validation of combination of control points, weights, and basis functions
        new{eltype(spline_dimensions), typeof(control_points), typeof(weights), typeof(eval)}(spline_dimensions, control_points, weights, eval)
    end
end

# TODO: Add indexing functionality such that indexing a spline_grid implies indexing spline_grid.eval

function SplineGrid(spline_dimensions::Vector{<:SplineDimension}, dim_out::Integer)::SplineGrid
    # The size of the point grid on which the spline is evaluated
    size_eval_grid = (length(spline_dimension.sample_points) for spline_dimension in spline_dimensions)
    # The size of the grid of control points
    size_cp_grid = get_n_basis_functions.(spline_dimensions)
    # The actual control points
    control_points = zeros(size_cp_grid..., dim_out)
    # Preallocated memory for grid evaluation of the spline
    eval = zeros(size_eval_grid..., dim_out)
    # NURBS are not supported yet
    weights = nothing
    SplineGrid(spline_dimensions, control_points, weights, eval)
end

function evaluate!(spline_grid::SplineGrid)::Nothing
    # TODO
    return nothing
end
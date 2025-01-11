"""
    DefaultControlPoints(control_points)

A thin wrapper around an array of control point values.
"""
struct DefaultControlPoints{
    Nin,
    Nout,
    Tv <: AbstractFloat,
    V <: AbstractArray{Tv}
} <: AbstractControlPoints{Nin, Nout, Tv}
    control_points::V
    function DefaultControlPoints(control_points)
        new{
            ndims(control_points) - 1,
            size(control_points)[end],
            eltype(control_points),
            typeof(control_points)
        }(
            control_points
        )
    end
end

function Base.show(
        io::IO,
        mime::MIME"text/plain",
        control_points::DefaultControlPoints{Nin, Nout, Tv}
) where {Nin, Nout, Tv}
    cp_grid_size = size(control_points)[1:(end - 1)]
    println(io,
        "DefaultControlPoints for grid of size $cp_grid_size in ℝ$(super(string(Nout))) ($Tv).")
end

Base.lastindex(control_points::DefaultControlPoints) = length(control_points)

function Base.getindex(control_points::DefaultControlPoints, inds...)
    view(control_points.control_points, inds...)
end

function Base.setindex!(control_points::DefaultControlPoints, val, inds...)
    control_points.control_points[inds...] = val
end

function Base.copyto!(control_points::DefaultControlPoints,
        vals::Union{AbstractArray, Base.Broadcast.Broadcasted}, inds...)
    copyto!(control_points.control_points, vals)
end

"""
    LocalRefinement(
        dim_refinement,
        refinement_matrix,
        refinement_indices,
        refinement_values)

The data for performing a local refinement.

## Fields

  - `dim_refinement`: The dimension along which will be refined
  - `refinement_matrix`: The matrix that performs the refinement
  - `refinement_indices`: The indices of the control points that will be overwritten after the refinement
  - `refinement_values`: The new values of the control points that will be rewritten after refinement
"""
struct LocalRefinement{
    Nin,
    Nout,
    Tv <: AbstractFloat,
    Ti <: Integer,
    R <: RefinementMatrix,
    I <: AbstractMatrix{Ti},
    D <: AbstractMatrix{Tv}
}
    dim_refinement::Int
    refinement_matrix::R
    refinement_indices::I
    refinement_values::D
    function LocalRefinement(
            dim_refinement,
            refinement_matrix,
            refinement_indices,
            refinement_values
    )
        new{
            size(refinement_indices)[2],
            size(refinement_values)[2],
            eltype(refinement_values),
            eltype(refinement_indices),
            typeof(refinement_matrix),
            typeof(refinement_indices),
            typeof(refinement_values)
        }(
            dim_refinement,
            refinement_matrix,
            refinement_indices,
            refinement_values
        )
    end
end

"""
    LocallyRefinedControlPoints(
        control_points_base,
        control_points_refined,
        local_refinements)

All data required to perform multiple refinement steps and overwriting control point values along the way,
yielding a Truncated Hierarchical Basis (THB) spline.

## Fields

  - `control_points_base`: The densely defined control points at the basis of the hierarchy
  - `control_points_refined`: The intermediate and final control point arrays after applying the local refinements
  - `local_refinements`: Data for local refinement for each step in the hierarchy
"""
struct LocallyRefinedControlPoints{
    Nin,
    Nout,
    Tv <: AbstractFloat,
    Ti <: Integer,
    V <: AbstractArray{Tv},
    L <: LocalRefinement{Nin, Nout, Tv, Ti}
} <: AbstractControlPoints{Nin, Nout, Tv}
    control_points_base::V
    control_points_refined::Vector{V}
    local_refinements::Vector{L}
    function LocallyRefinedControlPoints(
            control_points_base,
            control_points_refined,
            local_refinements::Vector{<:LocalRefinement{Nin, Nout, Tv, Ti}}
    ) where {Nin, Nout, Tv, Ti}
        new{
            ndims(control_points_base) - 1,
            size(control_points_base)[end],
            Tv,
            Ti,
            typeof(control_points_base),
            eltype(local_refinements)
        }(
            control_points_base,
            control_points_refined,
            local_refinements
        )
    end
end

# Helper function for getindex and setindex! for LocallyRefinedControlPoints
function get_control_point_view(control_points::LocallyRefinedControlPoints{
        Nin, Nout}) where {Nin, Nout}
    (; control_points_base, local_refinements) = control_points
    n_cp_base = prod(size(control_points_base)[1:(end - 1)])
    ApplyArray(
        vcat,
        reshape(control_points_base, n_cp_base, Nout),
        (lr.refinement_values for lr in local_refinements)...
    )
end

function Base.getindex(control_points::LocallyRefinedControlPoints, i_cp, i_dimout)
    view(get_control_point_view(control_points), i_cp, i_dimout)
end

function Base.setindex!(control_points::LocallyRefinedControlPoints, vals, i_cp, i_dimout)
    if length(vals) == 1
        get_control_point_view(control_points)[i_cp, i_dimout] = only(vals)
    else
        get_control_point_view(control_points)[i_cp, i_dimout] .= vals
    end
end

function Base.copyto!(control_points::LocallyRefinedControlPoints, vals)
    get_control_point_view(control_points) .= vals
end

function Base.show(
        io::IO,
        mime::MIME"text/plain",
        control_points::LocallyRefinedControlPoints{Nin, Nout, Tv}
) where {Nin, Nout, Tv}
    (; local_refinements) = control_points
    cp_grid_size = size(control_points)[1:(end - 1)]
    println(io,
        "LocallyRefinedControlPoints for final grid of size $cp_grid_size in ℝ$(super(string(Nout))) ($Tv). Local refinements:")
    data = zip(
        map(
        lr -> (
            lr.dim_refinement,
            lr.refinement_matrix.n,
            lr.refinement_matrix.m,
            size(lr.refinement_indices)[1]
        ),
        local_refinements
    )...)
    data = hcat(collect.(collect(data))...)
    header = ["input dim.", "# control p. before",
        "# control p. after", "# activated control p."]
    pretty_table(io, data; header)
end

Base.ndims(control_points::AbstractControlPoints) = ndims(obtain(control_points))
Base.ndims(::Type{C}) where {Nin, C <: AbstractControlPoints{Nin}} = Nin + 1
Base.size(control_points::AbstractControlPoints) = size(obtain(control_points))
Base.length(control_points::AbstractControlPoints) = length(obtain(control_points))
Base.eltype(::AbstractControlPoints{Nin, Nout, Tv}) where {Nin, Nout, Tv} = Tv
Base.vec(control_points::AbstractControlPoints) = vec(obtain(control_points))

"""
Get the number of control points within the control point object.
"""
function get_n_control_points(control_points::AbstractControlPoints{Nin}) where {Nin}
    prod(size(obtain(control_points))[1:Nin])
end

function get_n_control_points(control_points::LocallyRefinedControlPoints{Nin}) where {Nin}
    (; control_points_base, local_refinements) = control_points
    n_control_points = prod(size(control_points_base)[1:Nin])
    for local_refinement in local_refinements
        n_control_points += size(local_refinement.refinement_indices)[1]
    end
    n_control_points
end

function get_n_control_points(spline_grid::AbstractSplineGrid)
    get_n_control_points(spline_grid.control_points)
end

@kernel function local_refinement_kernel(
        control_points,
        @Const(refinement_indices),
        @Const(refinement_values)
)
    i = @index(Global, Linear)
    Nout = size(control_points)[end]
    refinement_indices_ = refinement_indices[i, :]

    for j in 1:Nout
        control_points[refinement_indices_..., j] = refinement_values[i, j]
    end
end

evaluate!(control_points::AbstractControlPoints) = nothing

"""
Evaluate the locally refined control points. For each local refinement, first
apply the refinement matrix and then overwrite the desired control point values.
"""
function evaluate!(control_points::LocallyRefinedControlPoints)
    (; control_points_base, control_points_refined, local_refinements) = control_points
    backend = get_backend(control_points_base)
    kernel! = local_refinement_kernel(backend)

    for (i, local_refinement) in enumerate(local_refinements)
        (; dim_refinement,
        refinement_matrix,
        refinement_indices,
        refinement_values) = local_refinement
        cp_prev = (i == 1) ? control_points_base : control_points_refined[i - 1]
        cp_new = control_points_refined[i]
        mult!(
            cp_new,
            refinement_matrix,
            cp_prev,
            dim_refinement
        )
        kernel!(
            cp_new,
            refinement_indices,
            refinement_values,
            ndrange = size(refinement_indices)[1]
        )
    end
end

obtain(control_points::AbstractArray) = control_points
obtain(control_points::DefaultControlPoints) = control_points.control_points

function obtain(control_points::LocallyRefinedControlPoints)
    if isempty(control_points.control_points_refined)
        control_points.control_points_base
    else
        last(control_points.control_points_refined)
    end
end

# Helper function for setting up or extending default local refinement
function default_local_refinement(
        spline_grid::AbstractSplineGrid{Nin, Nout, HasWeights, Tv, Ti},
        dim_refinement::Integer
) where {Nin, Nout, HasWeights, Tv, Ti}
    backend = get_backend(spline_grid.eval)
    spline_grid, refinement_matrix = refine(spline_grid, dim_refinement)

    refinement_indices = allocate(backend, Ti, 0, Nin)
    refinement_values = allocate(backend, Tv, 0, Nout)

    local_refinement = LocalRefinement(
        dim_refinement,
        refinement_matrix,
        refinement_indices,
        refinement_values
    )

    spline_grid, local_refinement
end

"""
Set up a default refinement matrix for each input dimension in order
(a default refinement matrix means bisecting each non-trivial knot span).
After this, data can be added to `spline_grid.control_points <: LocallyRefinedControlPoints` from
the output `spline_grid` to achieve local refinement by overwriting control point values
resulting from the refinement matrix multiplication.
"""
function setup_default_local_refinement(
        spline_grid::AbstractSplineGrid{Nin, Nout, false, Tv, Ti}
) where {Nin, Nout, Tv, Ti}
    control_points_base = obtain(spline_grid.control_points)

    # First dimension
    dim_refinement = 1
    spline_grid, local_refinement = default_local_refinement(
        spline_grid,
        dim_refinement
    )
    control_points_refined = [obtain(spline_grid.control_points)]
    local_refinements = [local_refinement]

    # Other dimensions
    for dim_refinement in 2:Nin
        spline_grid, local_refinement = default_local_refinement(
            spline_grid,
            dim_refinement
        )
        push!(control_points_refined, obtain(spline_grid.control_points))
        push!(local_refinements, local_refinement)
    end

    locally_refined_control_points = LocallyRefinedControlPoints(
        control_points_base,
        control_points_refined,
        local_refinements
    )

    @set spline_grid.control_points = locally_refined_control_points
end

"""
Extend default local refinement, that is: for a spline grid that already has locally refined control points,
add default local refinement for each input dimension in order.
"""
function extend_default_local_refinement(spline_grid::AbstractSplineGrid{Nin}) where {Nin}
    (; control_points) = spline_grid
    @assert control_points isa LocallyRefinedControlPoints
    (; local_refinements) = control_points

    for dim_refinement in 1:Nin
        spline_grid, local_refinement = default_local_refinement(
            spline_grid,
            dim_refinement
        )
        push!(local_refinements, local_refinement)
    end

    spline_grid
end

"""
Given a refinement step (by default the last refinement),
add new control points that overwrite the result from the refinement matrix.
These new values are chosen such that the spline geometry does not change.

## Inputs

  - `control_points`: The locally refined control points where new control points will be activated
  - `refinement_indices`: The indices in the control point grid at the `refinement_index` level which
    will be overwritten
  - `refinement_index`: The index of the refinement after which new control points will be activated
"""
function activate_local_refinement!(
        control_points::LocallyRefinedControlPoints{Nin, Nout, Tv, Ti},
        refinement_indices::AbstractMatrix{Ti};
        refinement_index::Integer = length(control_points.local_refinements)
)::Nothing where {Nin, Nout, Tv, Ti}
    (; local_refinements) = control_points
    control_points_refined = control_points.control_points_refined[refinement_index]
    local_refinement = local_refinements[refinement_index]

    @assert size(refinement_indices)[2]==Nin "Number of indices per control point must match the number of input dimensions."

    refinement_indices_new = unique(
        vcat(
            local_refinement.refinement_indices,
            refinement_indices
        ),
        dims = 1
    )

    n_refinements = size(local_refinement.refinement_indices)[1]
    n_refinements_new = size(refinement_indices_new)[1]

    # TODO: Do this more efficiently and GPU proof
    refinement_values_added = hcat([control_points_refined[
                                        refinement_indices_new[i, :]..., :]
                                    for i in (n_refinements + 1):n_refinements_new]...)'

    refinement_values_new = vcat(
        local_refinement.refinement_values,
        refinement_values_added
    )

    local_refinements[refinement_index] = LocalRefinement(
        local_refinement.dim_refinement,
        local_refinement.refinement_matrix,
        refinement_indices_new,
        refinement_values_new
    )

    return nothing
end

"""
Convenience wrapper of `activate_local_refinement!` for a spline grid
"""
function activate_local_refinement!(
        spline_grid::AbstractSplineGrid, args...; kwargs...)::Nothing
    activate_local_refinement!(
        spline_grid.control_points, args...; kwargs...)
end

"""
Instead of supplying explicit indices of the control points for local refinement activation,
supply a range of indices per dimension.
"""
function activate_local_control_point_range!(
        spline_grid::AbstractSplineGrid{Nin, Nout, HasWeights, Tv, Ti},
        ranges::Vararg{UnitRange}
)::Nothing where {Nin, Nout, HasWeights, Tv, Ti}
    @assert length(ranges) == Nin
    local_refinement_indices = [Ti.(collect(x)) for x in Iterators.product(ranges...)]
    local_refinement_indices = reduce(vcat, local_refinement_indices')
    activate_local_refinement!(spline_grid, local_refinement_indices)
end
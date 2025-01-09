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

Base.ndims(::Type{C}) where {Nin, C <: DefaultControlPoints{Nin}} = Nin + 1
Base.lastindex(control_points::DefaultControlPoints) = length(control_points)

function Base.getindex(control_points::DefaultControlPoints, inds...)
    view(control_points.control_points, inds...)
end

function Base.setindex!(control_points::DefaultControlPoints, val, inds...)
    control_points.control_points[inds...] = val
end

function Base.copyto!(control_points::DefaultControlPoints, vals, inds...)
    copyto!(control_points.control_points, vals)
end

function set_control_points!(
        control_points::DefaultControlPoints,
        values::AbstractArray
)
    copyto!(control_points.control_points, values)
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
    header = ["input dimension", "# control points before",
        "# control points after", "# overwritten control points"]
    pretty_table(io, data; header)
end

Base.ndims(control_points::AbstractControlPoints) = ndims(obtain(control_points))
Base.size(control_points::AbstractControlPoints) = size(obtain(control_points))
Base.length(control_points::AbstractControlPoints) = length(obtain(control_points))
Base.eltype(::AbstractControlPoints{Nin, Nout, Tv}) where {Nin, Nout, Tv} = Tv
Base.vec(control_points::AbstractControlPoints) = vec(obtain(control_points))

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

function set_control_points!(
        control_points::LocallyRefinedControlPoints{Nin, Nout},
        values::AbstractMatrix
)::Nothing where {Nin, Nout}
    (; control_points_base, local_refinements) = control_points
    n_control_points = get_n_control_points(control_points)
    n_input, n_output = size(values)
    @assert n_input==n_control_points "There are $n_control_points control points in this object, got $n_input."
    @assert n_output==Nout "These control points live in ℝ$(super(string(Nout))), got control points living in ℝ$(super(string(N_output)))."

    n_control_points_base = prod(size(control_points_base)[1:Nin])
    control_points_base .= reshape(
        view(values, 1:n_control_points_base, :), size(control_points_base))

    i_start = n_control_points_base + 1

    for local_refinement in local_refinements
        n_control_points_refinement = size(local_refinement.refinement_indices)[1]
        local_refinement.refinement_values .= view(
            values, i_start:(i_start + n_control_points_refinement - 1), :)
        i_start += n_control_points_refinement
    end
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
function KernelAbstractions.get_backend(control_points::AbstractControlPoints)
    get_backend(obtain(control_points))
end

function obtain(control_points::LocallyRefinedControlPoints)
    if isempty(control_points.control_points_refined)
        control_points.control_points_base
    else
        last(control_points.control_points_refined)
    end
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
    V = typeof(control_points_base)
    control_points_refined = V[]
    backend = get_backend(control_points_base)

    local local_refinements
    for dim in 1:Nin
        spline_grid, refinement_matrix = refine(spline_grid, dim)
        push!(control_points_refined, obtain(spline_grid.control_points))

        refinement_indices = allocate(backend, Ti, 0, Nin)
        refinement_values = allocate(backend, Tv, 0, Nout)

        local_refinement = LocalRefinement(
            dim,
            refinement_matrix,
            refinement_indices,
            refinement_values
        )

        if dim == 1
            local_refinements = [local_refinement]
        else
            push!(local_refinements, local_refinement)
        end
    end

    locally_refined_control_points = LocallyRefinedControlPoints(
        control_points_base,
        control_points_refined,
        local_refinements
    )

    @set spline_grid.control_points = locally_refined_control_points
end

function activate_local_refinement!(
        control_points::LocallyRefinedControlPoints{Nin, Nout, Tv, Ti},
        refinement_indices::AbstractMatrix{Ti};
        refinement_index::Integer = length(control_points.local_refinements)
)::Nothing where {Nin, Nout, Tv, Ti}
    (; local_refinements) = control_points
    control_points_refined = control_points.control_points_refined[refinement_index]
    local_refinement = local_refinements[refinement_index]

    @assert size(refinement_indices)[2]==Nin "Number of indices must match the number of input dimensions."

    refinement_indices_new = unique(
        vcat(
            local_refinement.refinement_indices,
            refinement_indices
        ),
        dims = 1
    )

    n_refinements = size(local_refinement.refinement_indices)[1]
    n_refinements_new = size(refinement_indices_new)[1]

    n_control_points_new = n_refinements_new - n_refinements

    # TODO: Do this more efficiently and GPU proof
    refinement_values_added = vcat([control_points_refined[
                                        refinement_indices_new[i, :]..., :]
                                    for i in (n_refinements + 1):n_refinements_new]...)

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

function activate_local_refinement!(
        spline_grid::AbstractSplineGrid, args...; kwargs...)::Nothing
    activate_local_refinement!(
        spline_grid.control_points, args...; kwargs...)
end
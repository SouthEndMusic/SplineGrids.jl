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

Base.ndims(control_points::AbstractControlPoints) = ndims(obtain(control_points))
Base.size(control_points::AbstractControlPoints) = size(obtain(control_points))
Base.length(control_points::AbstractControlPoints) = length(obtain(control_points))
Base.eltype(::AbstractControlPoints{Nin, Nout, Tv}) where {Nin, Nout, Tv} = Tv
Base.vec(control_points::AbstractControlPoints) = vec(obtain(control_points))

@kernel function local_refinement_kernel(
        control_points,
        @Const(refinement_indices),
        @Const(refinement_values)
)
    i = @index(Global, Linear)
    Nout = size(control_points)[end]
    refinement_indices_ = refinement_indices[i]

    for j in 1:Nout
        control_points[refinement_indices_..., j] = refinement_values[i, j]
    end
end

evaluate!(control_points::AbstractControlPoints) = nothing

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
            refinement_values
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

function refine_all_dimensions(
        spline_grid::AbstractSplineGrid{Nin, Nout, false, Tv, Ti}
) where {Nin, Nout, Tv, Ti}
    (; control_points) = spline_grid
    V = typeof(control_points)
    backend = get_backend(control_points)

    control_points_base = control_points
    control_points_refined = V[]

    local local_refinements
    for dim in 1:Nin
        spline_grid, refinement_matrix = refine(spline_grid, dim)
        push!(control_points_refined, spline_grid.control_points)

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
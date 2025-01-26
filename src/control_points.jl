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
        ::MIME"text/plain",
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
        dims_refinement,
        refinement_matrices,
        refinement_indices,
        refinement_values)

The data for performing a local refinement; one or more refinement matrices and the to be overwritten control point
values after multiplication with the refinement matrices along the specified dimensions.

## Fields

  - `dims_refinement`: The dimensions along which will be refined
  - `refinement_matrices`: The matrices that perform the refinement
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
    dims_refinement::Vector{Int}
    refinement_matrices::Vector{R}
    refinement_indices::I
    refinement_values::D
    function LocalRefinement(
            dims_refinement,
            refinement_matrices,
            refinement_indices,
            refinement_values
    )
        new{
            size(refinement_indices)[2],
            size(refinement_values)[2],
            eltype(refinement_values),
            eltype(refinement_indices),
            eltype(refinement_matrices),
            typeof(refinement_indices),
            typeof(refinement_values)
        }(
            dims_refinement,
            refinement_matrices,
            refinement_indices,
            refinement_values
        )
    end
end

function Adapt.adapt(backend::Backend, local_refinement::LocalRefinement)
    if get_backend(local_refinement.refinement_indices) == backend
        local_refinement
    else
        LocalRefinement(
            local_refinement.dims_refinement,
            map(refmat -> adapt(backend, refmat), local_refinement.refinement_matrices),
            adapt(backend, local_refinement.refinement_indices),
            adapt(backend, local_refinement.refinement_values)
        )
    end
end

# Set up a LocalRefinement for the base control points
function LocalRefinement(
        control_points_base::DefaultControlPoints{Nin, Nout, Tv},
        ::LocalRefinement{Nin, Nout, Tv, Ti, R, I, D}
)::LocalRefinement{Nin, Nout, Tv, Ti, R, I, D} where {Nin, Nout, Tv, Ti, R, I, D}
    backend = get_backend(control_points_base)
    n_control_points = get_n_control_points(control_points_base)
    refinement_indices = collect_indices(
        adapt(backend,
            reshape(
                collect(CartesianIndices(size(control_points_base)[1:(end - 1)])),
                n_control_points
            )
        ),
        Ti
    )
    refinement_values = copy(
        reshape(
        obtain(control_points_base),
        n_control_points,
        Nout
    )
    )
    LocalRefinement(
        Int[],
        R[],
        refinement_indices,
        refinement_values
    )
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
    control_points_refined::Vector{V}
    local_refinements::Vector{L}
    function LocallyRefinedControlPoints(
            control_points_refined,
            local_refinements::Vector{<:LocalRefinement{Nin, Nout, Tv, Ti}}
    ) where {Nin, Nout, Tv, Ti}
        control_points_base = first(control_points_refined)
        new{
            ndims(control_points_base) - 1,
            size(control_points_base)[end],
            Tv,
            Ti,
            typeof(control_points_base),
            eltype(local_refinements)
        }(
            control_points_refined,
            local_refinements
        )
    end
end

# Helper function for getindex and setindex! for LocallyRefinedControlPoints
function get_control_point_view(control_points::LocallyRefinedControlPoints{
        Nin, Nout}) where {Nin, Nout}
    (; local_refinements) = control_points
    ApplyArray(
        vcat,
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

function Base.copyto!(control_points::LocallyRefinedControlPoints,
        vals::Union{AbstractArray, Base.Broadcast.Broadcasted})
    get_control_point_view(control_points) .= vals
end

function Base.show(
        io::IO,
        ::MIME"text/plain",
        control_points::LocallyRefinedControlPoints{Nin, Nout, Tv}
) where {Nin, Nout, Tv}
    (; local_refinements) = control_points
    cp_grid_size = size(control_points)[1:(end - 1)]
    println(io,
        "LocallyRefinedControlPoints for final grid of size $cp_grid_size in ℝ$(super(string(Nout))) ($Tv). Local refinements:")

    local_refinement_base = first(local_refinements)

    input_dims = Number[NaN]
    n_cp_before = Number[NaN]
    n_cp_after = Number[NaN]
    n_cp_activated = Int[size(local_refinement_base.refinement_indices, 1)]

    for local_refinement in local_refinements[2:end]
        (; dims_refinement, refinement_matrices, refinement_indices) = local_refinement
        for (i, (input_dim, ref_mat)) in enumerate(zip(
            dims_refinement, refinement_matrices))
            push!(input_dims, input_dim)
            push!(n_cp_before, ref_mat.n)
            push!(n_cp_after, ref_mat.m)
            push!(n_cp_activated,
                (i == length(dims_refinement)) ? size(refinement_indices, 1) : 0)
        end
    end

    data = hcat(input_dims, n_cp_before, n_cp_after, n_cp_activated)
    header = ["input dim.", "# c.p. before", "# c.p. after", "# activated c.p."]
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
    (; local_refinements) = control_points
    n_control_points = 0
    for local_refinement in local_refinements
        n_control_points += size(local_refinement.refinement_indices, 1)
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

    Nin = ndims(control_points) - 1
    Nout = size(control_points)[end]

    indices = ntuple(dim_in -> refinement_indices[i, dim_in], Nin)

    for dim_out in 1:Nout
        control_points[indices..., dim_out] = refinement_values[i, dim_out]
    end
end

evaluate!(control_points::AbstractControlPoints) = nothing

"""
Evaluate the locally refined control points. For each local refinement, first
apply the refinement matrix and then overwrite the desired control point values.
"""
function evaluate!(control_points::LocallyRefinedControlPoints)
    (; control_points_refined, local_refinements) = control_points
    backend = get_backend(first(control_points_refined))
    kernel! = local_refinement_kernel(backend)

    for (i, local_refinement) in enumerate(local_refinements)
        (; dims_refinement,
        refinement_matrices,
        refinement_indices,
        refinement_values) = local_refinement
        cp_new = control_points_refined[i]
        if i > 1
            cp_prev = control_points_refined[i - 1]
            mult!(
                cp_new,
                Tuple(refinement_matrices),
                cp_prev,
                Tuple(dims_refinement)
            )
        end
        if !isempty(refinement_indices)
            kernel!(
                cp_new,
                refinement_indices,
                refinement_values,
                ndrange = (size(refinement_indices, 1),)
            )
            synchronize(backend)
        end
    end
end

obtain(control_points::AbstractArray) = control_points
obtain(control_points::DefaultControlPoints) = control_points.control_points

function obtain(control_points::LocallyRefinedControlPoints)
    last(control_points.control_points_refined)
end

"""
    add_default_local_refinement(spline_grid)

Refine in the default way in every dimension, i.e. bisecting every non-trivial knot span.
Yields a spline grid with a `LocallyRefinedControlPoints` object for the `control_points` field.
"""
function add_default_local_refinement(
        spline_grid::AbstractSplineGrid{
        Nin, Nout, HasWeights, Tv, Ti}
) where {Nin, Nout, HasWeights, Tv, Ti}
    (; spline_dimensions, control_points) = spline_grid
    backend = get_backend(control_points)

    # First dimension
    spline_dimension_new, refinement_matrix = refine(first(spline_dimensions))
    refinement_matrices = [refinement_matrix]
    spline_dimensions_new = ntuple(
        dim_in -> (dim_in == 1) ? spline_dimension_new : spline_dimensions[dim_in], Nin)

    # other dimensions
    for (dim_refinement, spline_dimension) in enumerate(spline_dimensions)
        (dim_refinement == 1) && continue
        spline_dimension_new, refinement_matrix = refine(spline_dimension)
        push!(refinement_matrices, refinement_matrix)
        spline_dimensions_new = ntuple(
            dim_in -> (dim_in == dim_refinement) ? spline_dimension_new :
                      spline_dimensions_new[dim_in],
            Nin)
    end

    evaluate!.(spline_dimensions_new)

    control_points_refined_new = KernelAbstractions.zeros(
        backend, Tv, get_control_point_grid_size(spline_dimensions_new)..., Nout)

    dims_refinement = Tuple(1:Nin)

    mult!(
        control_points_refined_new,
        Tuple(refinement_matrices),
        obtain(control_points),
        dims_refinement
    )

    local_refinement = LocalRefinement(
        collect(dims_refinement),
        refinement_matrices,
        allocate(backend, Ti, 0, Nin),
        allocate(backend, Tv, 0, Nout)
    )

    control_points_new = if control_points isa DefaultControlPoints
        local_refinement_base = LocalRefinement(control_points, local_refinement)
        LocallyRefinedControlPoints(
            [obtain(control_points), control_points_refined_new],
            [local_refinement_base, local_refinement]
        )
    else # control_points isa LocallyRefinedControlPoints
        push!(control_points.control_points_refined, control_points_refined_new)
        push!(control_points.local_refinements, local_refinement)
        control_points
    end

    setproperties(spline_grid;
        spline_dimensions = spline_dimensions_new,
        control_points = control_points_new
    )
end

@kernel function refinement_values_new_kernel(
        refinement_values_new,
        @Const(refinement_values),
        @Const(control_points_refined),
        @Const(refinement_indices_new)
)
    i = @index(Global, Linear)

    Nin = size(refinement_indices_new, 2)
    Nout = size(refinement_values_new, 2)
    n_refinement_values = size(refinement_values, 1)

    if i ≤ n_refinement_values
        for dim_out in 1:Nout
            refinement_values_new[i, dim_out] = refinement_values[i, dim_out]
        end
    else
        cp_grid_size = prod(size(control_points_refined)[1:(end - 1)])
        lin_cp_index_base = refinement_indices_new[i, 1]
        size_prod = 1
        for dim_in in 2:Nin
            size_prod *= size(control_points_refined, dim_in - 1)
            lin_cp_index_base += (refinement_indices_new[i, dim_in] - 1) * size_prod
        end
        for dim_out in 1:Nout
            lin_cp_idx = lin_cp_index_base + (dim_out - 1) * cp_grid_size
            refinement_values_new[i, dim_out] = control_points_refined[lin_cp_idx]
        end
    end
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
    backend = get_backend(control_points)

    @assert size(refinement_indices)[2]==Nin "Number of indices per control point must match the number of input dimensions."

    # Find unique refinement indices on CPU (https://github.com/JuliaGPU/CUDA.jl/issues/991)
    refinement_indices_new = adapt(
        backend,
        unique(
            adapt(
                CPU(),
                vcat(
                    local_refinement.refinement_indices,
                    refinement_indices
                )),
            dims = 1
        )
    )

    n_refinement_values_new = size(refinement_indices_new, 1)

    refinement_values_new = allocate(backend, Tv, n_refinement_values_new, Nout)
    refinement_values_new_kernel(backend)(
        refinement_values_new,
        local_refinement.refinement_values,
        control_points_refined,
        refinement_indices_new,
        ndrange = (n_refinement_values_new,)
    )
    synchronize(backend)

    local_refinements[refinement_index] = LocalRefinement(
        local_refinement.dims_refinement,
        local_refinement.refinement_matrices,
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

"""
    error_informed_local_refinement!(
        spline_grid::AbstractSplineGrid{Nin, Nout, HasWeights, Tv, Ti},
        error::AbstractArray;
        threshold::Union{Number, Nothing} = nothing
        )::Nothing where {Nin, Nout, HasWeights, Tv, Ti}

Refine the last level of the locally refined spline grid informed by the `error` array which has the same
shape as `spline_grid.eval`. This is done by:

  - mapping the error back onto the control points by using the adjoint of the refinement matrices multiplication
  - summing over the output dimensions to obtain a single number per control point stored in `control_grid_error`
  - activating each control point whose value is bigger than `threshold`

`threshold` can be explicitly provided but by default it is given by the mean of `control_grid_error`.
"""
function error_informed_local_refinement!(
        spline_grid::AbstractSplineGrid{Nin, Nout, HasWeights, Tv, Ti},
        error::AbstractArray;
        threshold::Union{Number, Nothing} = nothing
)::Nothing where {Nin, Nout, HasWeights, Tv, Ti}
    @assert size(error)==size(spline_grid.eval) "The error array must have the same size as the eval array."

    # Map the error array back onto the control point array
    control_points_error = zero(obtain(spline_grid.control_points))
    evaluate_adjoint!(spline_grid; eval = error, control_points = control_points_error)

    # Sum over output dimensions
    control_grid_error = dropdims(sum(control_points_error, dims = Nin + 1), dims = Nin + 1)

    # Default threshold if not provided
    if isnothing(threshold)
        threshold = sum(control_grid_error) / length(control_grid_error)
    end

    # Deduce the refinement indices
    refinement_indices_ci = findall(>(threshold), control_grid_error)
    refinement_indices = collect_indices(refinement_indices_ci, Ti)
    activate_local_refinement!(spline_grid, refinement_indices)

    return nothing
end

"""
Perform `deactivate_control_points` for every refinement level except the last one in reverse order.
For more details see `deactivate_overwritten_control_points!(::LocallyRefinedControlPoints, ::Integer)`.
"""
function deactivate_overwritten_control_points!(control_points::LocallyRefinedControlPoints)
    for level in (length(control_points.local_refinements) - 1):-1:1
        deactivate_overwritten_control_points!(control_points, level)
    end
end

"""
    deactivate_overwritten_control_points!(
        control_points::LocallyRefinedControlPoints,
        local_refinement_level::Integer)::Nothing

Deactivate control points whose effect is completely overwritten. The procedure works as follows:

The 'forward' computation to process local refinement is as follows:

`B = (O₂ ∘ L ∘ O₁)(A)`,

where:

  - `A` is the control point grid at `local_refinement_level`
  - `B` is the control point grid at `local_refinement_level + 1`
  - `O₁` is the overwriting operation of the active control points at `local_refinement_level`
  - `L` is the linear operation consisting of multiplying by refinement matrices along specified dimensions
  - `O₂` is the overwriting operation of the active control points at `local_refinement_level + 1`

To figure out whether for any overwriting element of `O₁` its effect is still present in the final array, we
perform the following computation:

`A = (L* ∘ O₂)(B)`,

where:

  - `B` is initialized with a special number `Flag`, which is a simple wrapper of a Boolean effectively saying
    whether this number is relevant or not. `B` is initialized with all `Flag(true)`
  - The overwriting values of `O₂` are initialized with all `Flag(false)`
  - `L*` is the adjoint version of `L`

After this computation we can read of in `A` at the overwriting locations of `O₁` whether each number is still having an
effect on `B` or not.
"""
function deactivate_overwritten_control_points!(
        control_points::LocallyRefinedControlPoints,
        local_refinement_level::Integer
)::Nothing
    (; local_refinements, control_points_refined) = control_points
    @assert 1 ≤ local_refinement_level ≤ length(local_refinements) - 1

    local_refinement = local_refinements[local_refinement_level]
    local_refinement_next = local_refinements[local_refinement_level + 1]

    backend = get_backend(control_points)

    n_refinement_next = size(local_refinement_next.refinement_indices, 1)
    refinement_values_next = KernelAbstractions.zeros(backend, Flag, n_refinement_next)
    control_points_next = KernelAbstractions.ones(
        backend, Flag, single_output_size(control_points_refined[local_refinement_level + 1]))

    local_refinement_kernel(backend)(
        control_points_next,
        local_refinement_next.refinement_indices,
        refinement_values_next,
        ndrange = (n_refinement_next,)
    )
    synchronize(backend)

    control_points_ = KernelAbstractions.zeros(
        backend, Flag, single_output_size(control_points_refined[local_refinement_level]))

    mult_adjoint!(
        control_points_,
        Tuple(local_refinement_next.refinement_matrices),
        control_points_next,
        Tuple(local_refinement_next.dims_refinement)
    )

    n_refinement_ = size(local_refinement.refinement_indices, 1)
    refinement_values_ = KernelAbstractions.zeros(backend, Flag, n_refinement_)

    local_refinement_adjoint_kernel(backend)(
        refinement_values_,
        control_points_,
        local_refinement.refinement_indices,
        ndrange = (n_refinement_,)
    )
    synchronize(backend)

    where_not_overwritten = findall(f -> f.flag, refinement_values_)

    refinement_indices_new = local_refinement.refinement_indices[where_not_overwritten, :]
    refinement_values_new = local_refinement.refinement_values[where_not_overwritten, :]
    local_refinements[local_refinement_level] = setproperties(
        local_refinement;
        refinement_indices = refinement_indices_new,
        refinement_values = refinement_values_new
    )
    return nothing
end
function set_unit_cp_grid!(control_points::AbstractArray)::Nothing
    Nin = ndims(control_points) - 1
    Nout = size(control_points)[end]
    cp_grid_size = size(control_points)
    indices = ntuple(_ -> :, Nin)
    for n in 1:min(Nin, Nout)
        n_coordinates = cp_grid_size[n]
        coordinates = collect(range(0, 1, length = n_coordinates))
        coordinates = reshape(
            coordinates, ntuple(i -> i == n ? cp_grid_size[i] : 1, Nin))
        control_points[indices..., n] .= repeat(
            coordinates, inner = ntuple(i -> i == n ? 1 : cp_grid_size[i], Nin))
    end
    return nothing
end

# Get the index i of each sample point t in the knot vector such 
# such that t ∈ [knots_all[i], knots_all[i + 1])
function set_sample_indices!(spline_dimension::AbstractSplineDimension)::Nothing
    backend = get_backend(spline_dimension)
    set_sample_indices_kernel(backend)(
        spline_dimension.sample_indices,
        spline_dimension.sample_points,
        spline_dimension.knot_vector.knots_all,
        spline_dimension.degree,
        ndrange = size(spline_dimension.sample_indices)
    )
    synchronize(backend)
end

# The number of basis functions in this spline dimension
function get_n_basis_functions(spline_dimension::AbstractSplineDimension)
    length(spline_dimension.knot_vector.knots_all) - spline_dimension.degree - 1
end

# Get the size of the block of control points that each output of the spline
# depends on
function get_cp_kernel_size(spline_dimensions::NTuple{
        Nin, <:AbstractSplineDimension})::NTuple{
        Nin, Int} where {Nin}
    ntuple(dim_in -> spline_dimensions[dim_in].degree + 1, Nin)
end

# The size of the grid on the domain of the spline where the 
# spline is evaluated
function get_sample_grid_size(spline_dimensions::NTuple{
        Nin, <:AbstractSplineDimension})::NTuple{
        Nin, Int} where {Nin}
    ntuple(n -> length(spline_dimensions[n].sample_points), Nin)
end

function get_control_point_grid_size(spline_dimensions::NTuple{
        Nin, <:AbstractSplineDimension})::NTuple{
        Nin, Int} where {Nin}
    get_n_basis_functions.(spline_dimensions)
end

function insert(v::AbstractVector{T}, i::Integer, x) where {T}
    backend = get_backend(v)
    out = allocate(backend, T, length(v) + 1)
    insert_kernel(backend)(out, v, i, x, ndrange = size(out))
    synchronize(backend)
    out
end

function shape_name(Nin::Integer)::String
    Nin == 1 && return "curve"
    Nin == 2 && return "surface"
    Nin == 3 && return "volume"
    return "hyper ($Nin) volume"
end

base_name(::AbstractSplineGrid) = "SplineGrid"
base_name(::AbstractNURBSGrid) = "NURBSGrid"

function Base.show(
        io::IO,
        mime::MIME"text/plain",
        spline_grid::AbstractSplineGrid{Nin, Nout, HasWeights, Tv}
) where {Nin, Nout, HasWeights, Tv}
    (; spline_dimensions, control_points) = spline_grid
    header = ["input dimension", "degree", "# basis functions", "# sample points"]
    data = zip(
        map(
        input -> begin
            (n, sd) = input
            (
                n,
                sd.degree,
                get_n_basis_functions(sd),
                length(sd.sample_points)
            )
        end,
        enumerate(spline_dimensions)
    )...)

    data = hcat(collect.(collect(data))...)
    intro = "$(base_name(spline_grid)) $(shape_name(Nin)) with outputs in ℝ$(super(string(Nout))) ($Tv)"

    println(io, intro, "\n", repeat("-", length(intro)))
    println(io, "* Properties per dimension:")
    pretty_table(io, data; header)
    println(io, "* Control points:")
    show(io, mime, control_points)
end

is_nurbs(::Any) = false
is_nurbs(::AbstractNURBSGrid) = true

function KernelAbstractions.get_backend(spline_dimension::AbstractSplineDimension)
    get_backend(spline_dimension.eval)
end

function KernelAbstractions.get_backend(control_points::AbstractControlPoints)
    get_backend(obtain(control_points))
end

function Adapt.adapt(
        backend::Backend,
        knot_vector::AbstractKnotVector{Ti, Tv}
)::AbstractKnotVector{Ti, Tv} where {Ti, Tv}
    if backend == get_backend(knot_vector.knot_values)
        knot_vector
    else
        KnotVector(
            adapt(backend, knot_vector.knot_values),
            adapt(backend, knot_vector.multiplicities)
        )
    end
end

function Adapt.adapt(
        backend::Backend,
        spline_dimension::AbstractSplineDimension{Ti, Tv}
)::AbstractSplineDimension{Ti, Tv} where {Ti, Tv}
    if backend == get_backend(spline_dimension)
        spline_dimension
    else
        SplineDimension(
            spline_dimension.degree,
            spline_dimension.max_derivative_order,
            adapt(backend, spline_dimension.knot_vector),
            adapt(backend, spline_dimension.sample_points),
            adapt(backend, spline_dimension.sample_indices),
            adapt(backend, spline_dimension.eval),
            adapt(backend, spline_dimension.eval_prev)
        )
    end
end

function Adapt.adapt(
        backend::Backend,
        control_points::AbstractControlPoints
)
    if get_backend(control_points) == backend
        control_points
    else
        if control_points isa DefaultControlPoints
            DefaultControlPoints(adapt(backend, control_points.control_points))
        else # control_points isa LocallyRefinedControlPoints
            LocallyRefinedControlPoints(
                map(cp -> adapt(backend, cp), control_points.control_points_refined),
                map(lr -> adapt(backend, lr), control_points.local_refinements)
            )
        end
    end
end

function Adapt.adapt(
        backend::Backend,
        spline_grid::AbstractSplineGrid{Nin, Nout, HasWeights, Tv, Ti}
)::AbstractSplineGrid{
        Nin, Nout, HasWeights, Tv, Ti} where {Nin, Nout, HasWeights, Tv, Ti}
    (; spline_dimensions) = spline_grid
    if backend == get_backend(first(spline_dimensions))
        spline_grid
    else
        SplineGrid(
            ntuple(i -> adapt(backend, spline_dimensions[i]), length(spline_dimensions)),
            adapt(backend, spline_grid.control_points),
            adapt(backend, spline_grid.weights),
            adapt(backend, spline_grid.eval)
        )
    end
end

function collect_indices(
        I::AbstractVector{<:CartesianIndex{Nin}},
        int_type::Type{Ti}
) where {Nin, Ti <: Integer}
    backend = get_backend(I)
    indices = KernelAbstractions.zeros(backend, int_type, length(I), Nin)
    collect_indices_kernel(backend)(indices, I, ndrange = size(I))
    synchronize(backend)
    indices
end

function single_output_size(array::AbstractArray)
    size_ = size(array)
    ndims_ = ndims(array)
    ntuple(dim -> (dim == ndims_) ? 1 : size_[dim], ndims_)
end

function get_row_extends(
        I,
        refmat_index_all,
        row_pointer_all,
        column_start_all,
        nzval_all
)
    Ndims = length(I)

    Ti = eltype(first(row_pointer_all))
    column_start = MVector{Ndims, Ti}(undef)
    n_columns = MVector{Ndims, Ti}(undef)

    for dim in 1:Ndims
        refmat_index = refmat_index_all[dim]
        if iszero(refmat_index)
            column_start[dim] = I[dim]
            n_columns[dim] = 1
        else
            column_start_, column_end = get_column_range(
                row_pointer_all[refmat_index],
                column_start_all[refmat_index],
                length(nzval_all[refmat_index]),
                I[dim]
            )
            column_start[dim] = column_start_
            n_columns[dim] = column_end - column_start_ + 1
        end
    end

    column_start, n_columns
end

# Flag is used as a simple number type to track which control points are completely overwritten
# in deactivate_overwritten_control_points!
struct Flag
    flag::Bool
end

Base.zero(::Type{Flag}) = Flag(false)
Base.one(::Type{Flag}) = Flag(true)
Base.convert(::Type{Flag}, x::Number) = Flag(Bool(x))

Base.:*(a::Flag, ::Number) = a
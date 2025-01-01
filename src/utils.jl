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
    set_sample_indices_kernel(get_backend(spline_dimension))(
        spline_dimension.sample_indices,
        spline_dimension.sample_points,
        spline_dimension.knot_vector.knots_all,
        spline_dimension.degree,
        ndrange = size(spline_dimension.sample_indices)
    )
    synchronize(backend)
end

function set_global_sample_indices!(
        sample_indices::AbstractArray{<:Integer},
        spline_dimensions::NTuple{
            Nin, <:AbstractSplineDimension},
        Nout::Integer)::Nothing where {Nin}
    cp_array_size = (get_control_point_grid_size(spline_dimensions)..., Nout)
    cp_indices = zeros(Int, Nin + 1)
    sample_indices_dims = [adapt(CPU(), spline_dim.sample_indices)
                           for spline_dim in spline_dimensions]
    for J in CartesianIndices(sample_indices)
        for dim in 1:Nin
            spline_dim = spline_dimensions[dim]
            cp_indices[dim] = sample_indices_dims[dim][J[dim]] -
                              spline_dim.degree - 1
        end
        sample_indices[J] = get_linear_index(cp_array_size, cp_indices)
    end
    return nothing
end

# Linear indices for control points per global sample point
# (Computed on CPU)
function get_global_sample_indices(
        spline_dimensions::NTuple{
            <:Any, <:AbstractSplineDimension},
        Nout::Integer)
    backend = get_backend(first(spline_dimensions))
    # The size of the point grid on which the spline is evaluated
    size_eval_grid = get_sample_grid_size(spline_dimensions)
    # Assumptions: dim_out = 0, I = (0,...,0)
    sample_indices = KernelAbstractions.zeros(CPU(), Int, size_eval_grid...)
    set_global_sample_indices!(sample_indices, spline_dimensions, Nout)
    adapt(backend, sample_indices)
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
    ntuple(n -> spline_dimensions[n].degree + 1, Nin)
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

# Outer product of n vectors, with thanks to Michael Abbott
function outer!(A::AbstractArray{T, N}, vs::Vararg{AbstractArray, N}) where {T, N}
    vecs = ntuple(n -> reshape(vs[n], ntuple(Returns(1), n - 1)..., :), N)
    broadcast!(*, A, vecs...)
end

# Type stable calculation of linear index
function get_linear_index(array_shape, indices)
    linear_idx = 1
    offset_local = 1

    for i in eachindex(indices)
        @inbounds linear_idx += (indices[i] - 1) * offset_local
        @inbounds offset_local *= array_shape[i]
    end

    return linear_idx
end

# Type stable calculate of the offset of the linear index for a
# location in the control point kernel given by `indices`
function get_offset(array_shape, indices)

    # Calculate flat index using column-major order (Julia's default)
    linear_idx = 0
    offset_local = 1

    for i in eachindex(indices)
        @inbounds linear_idx += indices[i] * offset_local
        @inbounds offset_local *= array_shape[i]
    end

    return linear_idx
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
        io::IO, mime::MIME"text/plain",
        spline_grid::AbstractSplineGrid{Nin, Nout}) where {Nin, Nout}
    (; spline_dimensions) = spline_grid
    header = ["input_dimension", "degree", "# basis functions", "# sample points"]
    data = zip(map(
        n -> (n, spline_dimensions[n].degree,
            SplineGrids.get_n_basis_functions(spline_dimensions[n]),
            length(spline_dimensions[n].sample_points)),
        eachindex(spline_dimensions))...)
    data = hcat(collect.(collect(data))...)

    println(
        io, "$(base_name(spline_grid)) $(shape_name(Nin)) with outputs in ℝ$(super(string(Nout))) with the following properties per dimension:")
    pretty_table(io, data; header)
end

is_nurbs(::Any) = false
is_nurbs(::AbstractNURBSGrid) = true

function KernelAbstractions.get_backend(spline_dimension::AbstractSplineDimension)
    get_backend(spline_dimension.eval)
end

function Adapt.adapt(backend::Backend,
        knot_vector::AbstractKnotVector{T})::AbstractKnotVector{T} where {T}
    KnotVector(
        adapt(backend, knot_vector.knot_values),
        adapt(backend, knot_vector.multiplicities)
    )
end

function Adapt.adapt(backend::Backend,
        spline_dimension::AbstractSplineDimension{T})::AbstractSplineDimension{T} where {T}
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

function Adapt.adapt(backend::Backend,
        spline_grid::AbstractSplineGrid{Nin, Nout, HasWeights, T}
)::AbstractSplineGrid{
        Nin, Nout, HasWeights, T} where {Nin, Nout, HasWeights, T}
    (; spline_dimensions) = spline_grid
    SplineGrid(
        ntuple(i -> adapt(backend, spline_dimensions[i]), length(spline_dimensions)),
        adapt(backend, spline_grid.control_points),
        adapt(backend, spline_grid.eval),
        adapt(backend, spline_grid.sample_indices),
        adapt(backend, spline_grid.basis_function_products);
        denominator = adapt(backend, spline_grid.denominator),
        weights = adapt(backend, spline_grid.weights)
    )
end
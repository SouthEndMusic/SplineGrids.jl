function set_unit_cp_grid!(control_points::AbstractArray)::Nothing
    N_in = ndims(control_points) - 1
    N_out = size(control_points)[end]
    cp_grid_size = size(control_points)
    indices = ntuple(_ -> :, N_in)
    for n in 1:min(N_in, N_out)
        n_coordinates = cp_grid_size[n]
        coordinates = collect(range(0, 1, length = n_coordinates))
        coordinates = reshape(
            coordinates, ntuple(i -> i == n ? cp_grid_size[i] : 1, N_in))
        control_points[indices..., n] .= repeat(
            coordinates, inner = ntuple(i -> i == n ? 1 : cp_grid_size[i], N_in))
    end
    return nothing
end

function set_global_sample_indices!(
        sample_indices::AbstractArray{<:Integer},
        spline_dimensions::NTuple{
            Nin, <:AbstractSplineDimension},
        control_points::AbstractArray)::Nothing where {Nin}
    cp_indices = zeros(Int, Nin + 1)
    for J in CartesianIndices(sample_indices)
        for dim in 1:Nin
            spline_dim = spline_dimensions[dim]
            cp_indices[dim] = spline_dim.sample_indices[J[dim]] -
                              spline_dim.degree - 1
        end
        sample_indices[J] = get_linear_index(size(control_points), cp_indices)
    end
    return nothing
end

# Linear indices for control points per global sample point
function get_global_sample_indices(
        spline_dimensions::NTuple{
            Nin, <:AbstractSplineDimension},
        control_points::AbstractArray) where {Nin}
    # The size of the point grid on which the spline is evaluated
    size_eval_grid = ntuple(n -> length(spline_dimensions[n].sample_points), Nin)
    # Assumptions: dim_out = 0, I = (0,...,0)
    sample_indices = zeros(Int, size_eval_grid...)
    set_global_sample_indices!(sample_indices, spline_dimensions, control_points)
    sample_indices
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
function outer!(A::AbstractArray{T, N}, vs::Vararg{SubArray, N}) where {T, N}
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

function reset_sample_indices!(spline_dimension::AbstractSplineDimension)::Nothing
    (; knot_vector, degree, sample_points, sample_indices) = spline_dimension
    map!(sample_point -> get_index(
            knot_vector,
            sample_point,
            degree),
        sample_indices,
        sample_points
    )
    return nothing
end

is_nurbs(::Any) = false
is_nurbs(::AbstractNURBSGrid) = true
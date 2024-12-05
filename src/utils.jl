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

function shape_name(Nin::Integer)::String
    Nin == 1 && return "curve"
    Nin == 2 && return "surface"
    Nin == 3 && return "volume"
    return "hyper ($Nin) volume"
end

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
        io, "SplineGrid $(shape_name(Nin)) in ℝ$(super(string(Nout))) with the following properties per dimension:")
    pretty_table(io, data; header)
end
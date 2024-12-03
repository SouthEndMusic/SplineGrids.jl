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
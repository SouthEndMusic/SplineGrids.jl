module SplineGridsMakieExt
using Accessors
using Adapt
using KernelAbstractions
using Makie
using SplineGrids

function SplineGrids.plot_basis(spline_grid::SplineGrid{2})
    fig = Figure()
    ax = Axis3(fig[1, 1], azimuth = -π / 2, elevation = π / 2, perspectiveness = 0.5)
    plot_basis!(ax, spline_grid)
    fig
end

function SplineGrids.plot_basis!(
        ax::Makie.Axis3,
        spline_grid::SplineGrid{2, Nout, Tv}
)::Nothing where {Nout, Tv}
    (; control_points, spline_dimensions) = spline_grid
    spline_grid = adapt(CPU(), spline_grid)
    control_points_new = deepcopy(control_points)
    @reset spline_grid.control_points = control_points_new
    backend = get_backend(spline_grid.control_points)

    n_control_points = get_n_control_points(spline_grid)
    values = KernelAbstractions.zeros(backend, Tv, n_control_points, Nout)

    eval_view = view(spline_grid.eval, :, :, 1)
    basis_functions = zero(eval_view)
    for i in 1:n_control_points
        values .= 0
        values[i, 1] = 1
        set_control_points!(spline_grid, values)
        evaluate!(control_points_new)
        evaluate!(spline_grid)
        broadcast!(max, basis_functions, basis_functions, eval_view)
    end
    surface!(
        ax,
        spline_dimensions[1].sample_points,
        spline_dimensions[2].sample_points,
        basis_functions,
        colormap = :coolwarm,
        colorrange = (0, 1)
    )
    return nothing
end

end # module SplineGridsMakieExt
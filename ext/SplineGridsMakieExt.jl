module SplineGridsMakieExt
using ConstructionBase
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
    spline_grid = adapt(CPU(), spline_grid)
    (; control_points, spline_dimensions) = spline_grid
    control_points_new = deepcopy(control_points)
    spline_grid = setproperties(spline_grid; control_points = control_points_new)

    eval_view = view(spline_grid.eval, :, :, 1)
    basis_functions = zero(eval_view)
    for i in 1:get_n_control_points(spline_grid)
        control_points_new .= 0
        if control_points_new isa DefaultControlPoints
            control_points_new[CartesianIndices(size(control_points_new)[1:2])[i], 1] = 1
        else
            control_points_new[i, 1] = 1
        end
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

module SplineGridsMakieExt
using ConstructionBase
using Adapt
using KernelAbstractions
using Makie
using SplineGrids

function SplineGrids.plot_basis(spline_grid::SplineGrid{2}; kwargs...)
    fig = Figure()
    plot_basis!(fig, spline_grid; kwargs...)
    fig
end

function SplineGrids.plot_basis!(
        fig::Figure, spline_grid::SplineGrid{2}; i = 1, j = 1, kwargs...)
    extents = ntuple(i -> spline_grid.spline_dimensions[i].knot_vector.extent, 2)
    ratio = (extents[2][2] - extents[2][1]) / (extents[1][2] - extents[1][1])
    ax = Axis3(
        fig[i, j], azimuth = -π / 2, elevation = π / 2,
        perspectiveness = 0.5, aspect = (1, ratio, 1); kwargs...)
    plot_basis!(ax, spline_grid)
end

function SplineGrids.plot_basis!(ax::Makie.Axis3, spline_grid::SplineGrid{2})::Nothing
    spline_grid = adapt(CPU(), spline_grid)
    (; control_points) = spline_grid
    control_points_new = deepcopy(control_points)
    spline_grid = setproperties(spline_grid; control_points = control_points_new)

    plot_basis!(ax, spline_grid, control_points_new)
    return nothing
end

function SplineGrids.plot_basis!(
        ax::Makie.Axis3,
        spline_grid::SplineGrid{2},
        control_points_new::DefaultControlPoints
)
    (; spline_dimensions) = spline_grid
    eval_view = view(spline_grid.eval, :, :, 1)
    basis_functions = zero(eval_view)
    CI = CartesianIndices(size(control_points_new)[1:2])

    for i in 1:get_n_control_points(spline_grid)
        control_points_new .= 0
        control_points_new[CI[i], 1] = 1
        evaluate!(control_points_new)
        evaluate!(spline_grid)
        broadcast!(max, basis_functions, basis_functions, eval_view)
    end
    surface!(
        ax,
        spline_dimensions[1].sample_points,
        spline_dimensions[2].sample_points,
        basis_functions,
        color = cbrt.(basis_functions),
        colormap = [:black, first(Makie.wong_colors())],
        colorrange = (0, 1)
    )
end

function SplineGrids.plot_basis!(
        ax::Makie.Axis3,
        spline_grid::SplineGrid{2, Nout, Tv},
        control_points_new::LocallyRefinedControlPoints
) where {Nout, Tv}
    (; local_refinements) = control_points_new
    (; spline_dimensions) = spline_grid
    eval_view = view(spline_grid.eval, :, :, 1)
    control_points_new .= 0
    basis_functions = zero(eval_view)
    wong_colors = Makie.wong_colors()
    colors = fill(Makie.RGBA{Tv}(0, 0, 0, 1), size(eval_view))

    for (level, local_refinement) in enumerate(local_refinements)
        (; refinement_values) = local_refinement
        basis_functions_level = zero(eval_view)
        for i in axes(refinement_values, 1)
            refinement_values .= 0
            refinement_values[i, 1] = 1
            evaluate!(control_points_new)
            evaluate!(spline_grid)
            broadcast!(max, basis_functions_level, basis_functions_level, eval_view)
        end
        cmap_level = cgrad([:black, wong_colors[level]])
        for I in eachindex(basis_functions)
            value = basis_functions[I]
            value_level = basis_functions_level[I]
            if value_level > value
                basis_functions[I] = value_level
                colors[I] = get(cmap_level, ∛(value_level))
            end
        end
    end
    surface!(
        ax,
        spline_dimensions[1].sample_points,
        spline_dimensions[2].sample_points,
        basis_functions,
        color = colors
    )
end

end # module SplineGridsMakieExt

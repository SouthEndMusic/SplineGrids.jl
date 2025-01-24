module SplineGridsLinearMapsExt
using SplineGrids
using LinearMaps

function LinearMaps.LinearMap(
        spline_grid::SplineGrid{Nin, Nout, Tv};
        derivative_order::NTuple{Nin, <:Integer} = ntuple(_ -> 0, Nin)
) where {Nin, Nout, Tv}
    SplineGrids.validate_partial_derivatives(
        spline_grid, derivative_order
    )

    LinearMap{Tv}(
        # In place evaluation control points -> spline grid
        (evaluation_flat, control_points_flat) -> evaluate!(
            spline_grid;
            derivative_order,
            control_points = reshape(control_points_flat, size(spline_grid.control_points)),
            eval = reshape(evaluation_flat, size(spline_grid.eval))
        ),
        # In place evaluation spline grid -> control points (adjoint)
        (control_points_flat, evaluation_flat) -> begin
            evaluate_adjoint!(
                spline_grid;
                derivative_order,
                control_points = reshape(
                    control_points_flat, size(spline_grid.control_points)),
                eval = reshape(evaluation_flat, size(spline_grid.eval))
            )
            # TODO: This must be the adjoint of evaluate!(control_points::LocallyRefinedControlPoints)
            # if applicable
            # control_points_flat .= vec(SplineGrids.obtain(spline_grid.control_points))
        end,
        length(spline_grid.eval),
        get_n_control_points(spline_grid) * Nout
    )
end

end # Module SplineGridsLinearMapsExt

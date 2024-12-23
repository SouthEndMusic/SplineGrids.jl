module SplineGridsLinearMapsExt
using SplineGrids
using LinearMaps

function LinearMaps.LinearMap(spline_grid::SplineGrids.AbstractSplineGrid{Nin};
        derivative_order::NTuple{Nin, <:Integer} = ntuple(_ -> 0, Nin)) where {Nin}
    SplineGrids.validate_partial_derivatives(
        spline_grid, derivative_order
    )

    LinearMap(
        # In place evaluation control points -> spline grid
        (evaluation_flat, control_points_flat) -> evaluate!(
            spline_grid;
            derivative_order,
            control_points = reshape(control_points_flat, size(spline_grid.control_points)),
            eval = reshape(evaluation_flat, size(spline_grid.eval))
        ),
        # In place evaluation spline grid -> control points (adjoint)
        (control_points_flat, evaluation_flat) -> evaluate_adjoint!(
            spline_grid;
            derivative_order,
            control_points = reshape(control_points_flat, size(spline_grid.control_points)),
            eval = reshape(evaluation_flat, size(spline_grid.eval))
        ),
        length(spline_grid.eval),
        length(spline_grid.control_points)
    )
end

end # Module SplineGridsLinearMapsExt
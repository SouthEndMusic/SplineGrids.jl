module SplineGridsLinearMapsExt
using SplineGrids
using LinearMaps

function LinearMaps.LinearMap(
        spline_grid::SplineGrid{Nin, Nout, Tv};
        derivative_order::NTuple{Nin, <:Integer} = ntuple(_ -> 0, Nin)
) where {Nin, Nout, Tv}
    (; control_points) = spline_grid
    SplineGrids.validate_partial_derivatives(
        spline_grid, derivative_order
    )

    LinearMap{Tv}(
        # In place evaluation control points -> spline grid
        (evaluation_flat, control_points_flat) -> begin
            copyto!(control_points, control_points_flat)
            evaluate!(control_points)
            evaluate!(
                spline_grid;
                derivative_order,
                eval = reshape(evaluation_flat, size(spline_grid.eval))
            )
        end,
        # In place evaluation spline grid -> control points (adjoint)
        (control_points_flat, evaluation_flat) -> begin
            evaluate_adjoint!(
                spline_grid;
                derivative_order,
                eval = reshape(evaluation_flat, size(spline_grid.eval))
            )
            evaluate_adjoint!(control_points)
            copyto!(control_points_flat, control_points)
        end,
        length(spline_grid.eval),
        get_n_control_points(spline_grid) * Nout
    )
end

end # Module SplineGridsLinearMapsExt

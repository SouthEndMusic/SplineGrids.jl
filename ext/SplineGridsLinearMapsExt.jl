module SplineGridsLinearMapsExt
using SplineGrids
using LinearMaps

function LinearMaps.LinearMap(spline_grid::SplineGrid)
    LinearMap(
        (evaluation_flat, control_points_flat) -> evaluate!(
            spline_grid;
            control_points = reshape(control_points_flat, size(spline_grid.control_points)),
            eval = reshape(evaluation_flat, size(spline_grid.eval))
        ),
        length(spline_grid.eval),
        length(spline_grid.control_points)
    )
end

end # Module SplineGridsLinearMapsExt
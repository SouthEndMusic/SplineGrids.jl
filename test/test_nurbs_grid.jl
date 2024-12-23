using SplineGrids

@testset "Circle" begin
    using SplineGrids
    using Plots

    n_control_points = 9
    degree = 2
    n_sample_points = 500
    dim_out = 2

    knot_values = [0, π / 2, π, 3π / 2, 2π]
    multiplicities = [3, 2, 2, 2, 3]
    knot_vector = KnotVector(knot_values, multiplicities)

    spline_dimension = SplineDimension(
        n_control_points, degree, n_sample_points; knot_vector)
    nurbs_grid = NURBSGrid(spline_dimension, dim_out)
    nurbs_grid.weights[2:2:end] .= 1 / sqrt(2)
    nurbs_grid.control_points .= [1 0;
                                  1 1;
                                  0 1;
                                  -1 1;
                                  -1 0;
                                  -1 -1;
                                  0 -1;
                                  1 -1;
                                  1 0]
    evaluate!(nurbs_grid)
    points_on_circle = eachrow(nurbs_grid.eval)
    @test allunique(points_on_circle[2:end])
    @test !any(point -> !(norm(point) ≈ 1), points_on_circle)
end
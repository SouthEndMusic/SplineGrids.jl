using SplineGrids
using KernelAbstractions
using Adapt

if "--gpu_backend" ∈ ARGS
    backend = CUDABackend()
else
    backend = CPU()
end

@testset "Validation" begin
    n_control_points = 4
    degree = 2
    n_sample_points = 100
    dim_out = 1
    max_derivative_order = 1

    spline_dimension = SplineDimension(
        n_control_points, degree, n_sample_points; max_derivative_order, backend)
    nurbs_grid = NURBSGrid(spline_dimension, dim_out)
    @test_throws "Computing derivatives of NURBS is currently not supported." evaluate!(
        nurbs_grid; derivative_order = (1,))
end

@testset "Circle" begin
    n_control_points = 9
    degree = 2
    n_sample_points = 500
    dim_out = 2

    knot_values = adapt(backend, Float32[0, π / 2, π, 3π / 2, 2π])
    multiplicities = adapt(backend, Int32[3, 2, 2, 2, 3])
    knot_vector = KnotVector(knot_values, multiplicities)

    spline_dimension = SplineDimension(
        n_control_points, degree, n_sample_points; knot_vector)
    nurbs_grid = NURBSGrid(spline_dimension, dim_out)
    nurbs_grid.weights[2:2:end] .= 1 / sqrt(2)
    copyto!(nurbs_grid.control_points, [1 0;
                                        1 1;
                                        0 1;
                                        -1 1;
                                        -1 0;
                                        -1 -1;
                                        0 -1;
                                        1 -1;
                                        1 0])
    evaluate!(nurbs_grid)
    points_on_circle = eachrow(adapt(CPU(), nurbs_grid.eval))
    @test allunique(points_on_circle[2:end])
    @test all(point -> point[1]^2 + point[2]^2 ≈ 1, points_on_circle)
end
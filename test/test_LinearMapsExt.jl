using SplineGrids
using LinearMaps
using IterativeSolvers
using Random
using KernelAbstractions
using Adapt

if "--gpu_backend" ∈ ARGS
    backend = CUDABackend()
else
    backend = CPU()
end

@testset "Least squares fitting" begin
    Random.seed!(1)

    n_control_points = (12, 10)
    degree = (3, 4)
    n_sample_points = (26, 73)
    Nout = 3

    spline_dimensions = SplineDimension.(
        n_control_points, degree, n_sample_points; backend, float_type = Float64)
    spline_grid = SplineGrid(spline_dimensions, Nout)
    control_points_rand = adapt(backend, rand(size(spline_grid.control_points)...))
    copyto!(spline_grid.control_points, control_points_rand)
    evaluate!(spline_grid)

    spline_grid_map = LinearMap(spline_grid)
    fit = lsqr(spline_grid_map, copy(vec(spline_grid.eval)))
    @test fit≈vec(control_points_rand) rtol=1e-5
end

@testset "Locally refined least squares fitting" begin
    Random.seed!(1)

    n_control_points = (12, 10)
    degree = (3, 4)
    n_sample_points = (26, 73)
    Nout = 3

    spline_dimensions = SplineDimension.(
        n_control_points, degree, n_sample_points; backend, float_type = Float64)
    spline_grid = SplineGrid(spline_dimensions, Nout)
    spline_grid = add_default_local_refinement(spline_grid)

    activate_local_control_point_range!(spline_grid, 1:4, 1:6)
    activate_local_control_point_range!(spline_grid, 1:6, 1:2)
    activate_local_control_point_range!(spline_grid, 9:10, 7:10)
    deactivate_overwritten_control_points!(spline_grid.control_points)

    n_control_points = get_n_control_points(spline_grid)
    control_points_rand = adapt(backend, rand(n_control_points, Nout))
    copyto!(spline_grid.control_points, control_points_rand)
    evaluate!(spline_grid.control_points)
    evaluate!(spline_grid)

    eval_target = copy(spline_grid.eval)
    spline_grid_map = LinearMap(spline_grid)
    fit = lsqr(spline_grid_map, vec(eval_target))
    @test fit≈vec(control_points_rand) rtol=1e-5
end

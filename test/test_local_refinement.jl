using SplineGrids
using KernelAbstractions

if "--gpu_backend" ∈ ARGS
    backend = CUDABackend()
else
    backend = CPU()
end

@testset "Local refinement" begin
    using SplineGrids: obtain

    n_control_points = (6, 6)
    degree = (2, 2)
    n_sample_points = (500, 500)
    dim_out = 3

    spline_dimensions = SplineDimension.(n_control_points, degree, n_sample_points; backend)
    spline_grid = SplineGrid(spline_dimensions, dim_out)
    @test get_n_control_points(spline_grid) == 36
    @test isnothing(evaluate!(spline_grid.control_points))

    spline_grid = setup_default_local_refinement(spline_grid)
    (; control_points) = spline_grid

    @test control_points isa LocallyRefinedControlPoints
    @test occursin(
        "LocallyRefinedControlPoints for final grid of size (10, 10) in ℝ³ (Float32). Local refinements:",
        sprint(io -> show(io, MIME"text/plain"(), spline_grid))
    )
    @test get_n_control_points(spline_grid) == 36

    control_points_before = copy(obtain(control_points))

    activate_local_control_point_range!(spline_grid, 1:4, 1:6)
    activate_local_control_point_range!(spline_grid, 1:6, 1:2)
    activate_local_control_point_range!(spline_grid, 9:10, 7:10)
    evaluate!(control_points)

    @test get_n_control_points(spline_grid) == 72
    @test obtain(control_points) ≈ control_points_before

    spline_grid = extend_default_local_refinement(spline_grid)

    empty!(control_points.control_points_refined)
    @test obtain(control_points) === control_points.control_points_base
end
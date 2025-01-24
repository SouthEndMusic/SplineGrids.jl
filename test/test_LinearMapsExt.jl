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
    copyto!(spline_grid.control_points, rand(size(spline_grid.control_points)...))
    evaluate!(spline_grid)

    spline_grid_map = LinearMap(spline_grid)
    fit = lsqr(spline_grid_map, vec(spline_grid.eval))
    @test fit≈vec(spline_grid.control_points) rtol=1e-5
end

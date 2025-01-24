using SplineGrids
using KernelAbstractions

if "--gpu_backend" ∈ ARGS
    backend = CUDABackend()
else
    backend = CPU()
end

n_control_points = (5, 6)
degree = (3, 2)
n_sample_points = (7, 9)
Nout = 2

spline_dimensions = SplineDimension.(
    n_control_points,
    degree,
    n_sample_points;
    backend
)
spline_grid = SplineGrid(spline_dimensions, Nout)

@testset "Control point grid size" begin
    @test SplineGrids.get_control_point_grid_size(spline_dimensions) == n_control_points
end

@testset "Show" begin
    @test startswith(sprint(io -> show(io, MIME"text/plain"(), spline_grid)),
        "SplineGrid surface with outputs in ℝ² (Float32)\n")
end

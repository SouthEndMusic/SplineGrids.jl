using SplineGrids
using Random
using KernelAbstractions

if "--gpu_backend" ∈ ARGS
    backend = CUDABackend()
else
    backend = CPU()
end

Random.seed!(1)

n_control_points = (5, 8, 6)
degree = (4, 2, 3)
n_sample_points = (15, 20, 25)
Nout = 2

spline_dimensions = SplineDimension.(
    n_control_points,
    degree,
    n_sample_points;
    distribution = :random,
    backend
)
spline_grid_ = SplineGrid(spline_dimensions, Nout)
copyto!(spline_grid_.control_points, rand(n_control_points..., Nout))
evaluate!(spline_grid_)

@testset "Knot insertion" begin
    spline_grid = deepcopy(spline_grid_)
    @test spline_grid_.eval !== spline_grid.eval # Make sure that these are different arrays

    for dim_refinement in eachindex(n_control_points)
        spline_grid, refinement_matrix = insert_knot!(spline_grid, dim_refinement, 0.25)
        evaluate!(spline_grid)
        @test spline_grid_.eval ≈ spline_grid.eval # Knot insertion didn't change geometry
        size_cp_grid = size(spline_grid.control_points)
        @test size(refinement_matrix) == (
            size_cp_grid[dim_refinement], size_cp_grid[dim_refinement] - 1)
    end
end

@testset "Refinement" begin
    spline_grid = deepcopy(spline_grid_)
    @test spline_grid_.eval !== spline_grid.eval # Make sure that these are different arrays

    size_cp_grid_old = size(spline_grid_.control_points)

    for dim_refinement in eachindex(n_control_points)
        spline_grid, refinement_matrix = refine!(spline_grid, dim_refinement)
        evaluate!(spline_grid)
        @test spline_grid_.eval ≈ spline_grid.eval # Knot insertion didn't change geometry
        size_cp_grid = size(spline_grid.control_points)
        @test size(refinement_matrix) == (
            size_cp_grid[dim_refinement], size_cp_grid_old[dim_refinement])
    end
end
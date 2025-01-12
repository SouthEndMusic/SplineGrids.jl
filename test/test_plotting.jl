using SplineGrids
using KernelAbstractions

if "--gpu_backend" ∈ ARGS
    backend = CUDABackend()
else
    backend = CPU()
end

# Avoid conflict with Plots.backend
backend_ = backend

using Plots

n_control_points = (5, 6)
degree = (3, 2)
n_sample_points = (32, 55)
max_derivative_order = 1

spline_dimensions = SplineDimension.(
    n_control_points,
    degree,
    n_sample_points;
    max_derivative_order,
    backend = backend_
)

@testset "Plot spline dimension" begin
    @test try
        plot(spline_dimensions[1])
        true
    catch e
        println("Error: $e")
        Base.show_backtrace(stderr, catch_backtrace())
        false
    end

    @test try
        plot(spline_dimensions[1]; derivative_order = 1)
        true
    catch e
        println("Error: $e")
        Base.show_backtrace(stderr, catch_backtrace())
        false
    end
end

dimensionalities = [
    ((1, Nout) for Nout in 1:4)...,
    ((2, Nout) for Nout in 1:3)...
]

for (Nin, Nout) in dimensionalities
    @testset "Plot spline grid, Nin = $Nin, Nout = $Nout" begin
        spline_grid = SplineGrid(spline_dimensions[1:Nin], Nout)
        @test try
            plot(spline_grid)
            true
        catch e
            println("Error: $e")
            Base.show_backtrace(stderr, catch_backtrace())
            false
        end

        nurbs_grid = NURBSGrid(spline_dimensions[1:Nin], Nout)
        @test try
            plot(nurbs_grid)
            true
        catch e
            println("Error: $e")
            Base.show_backtrace(stderr, catch_backtrace())
            false
        end
    end
end

@testset "MakieExt" begin
    using CairoMakie
    spline_grid = SplineGrid(spline_dimensions[1:2], 3)

    @test try
        plot_basis(spline_grid)
        true
    catch e
        println("Error: $e")
        Base.show_backtrace(stderr, catch_backtrace())
        false
    end

    spline_grid = setup_default_local_refinement(spline_grid)
    @test try
        plot_basis(spline_grid)
        true
    catch e
        println("Error: $e")
        Base.show_backtrace(stderr, catch_backtrace())
        false
    end
end
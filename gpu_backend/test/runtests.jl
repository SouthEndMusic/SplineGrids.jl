using Test
using CUDA: CUDABackend

push!(ARGS, "--gpu_backend")

test_dir = normpath(@__DIR__, "../../test")

test_data = [("Utilities", "test_utils.jl")
             ("KnotVector", "test_knot_vector.jl")
             ("SplineDimension", "test_spline_dimension.jl")
             ("SplineGrid", "test_spline_grid.jl")
             ("Refinement", "test_refinement.jl")
             ("NURBSGrid", "test_nurbs_grid.jl")
             ("Plotting", "test_plot_rec.jl")
             ("LinearMapsExt", "test_LinearMapsExt.jl") # Not supported for CUDA backend?
             ("EnzymeExt", "test_EnzymeExt.jl")]

for (test_name, test_file) in test_data
    @testset "$test_name" begin
        include(normpath(test_dir, test_file))
        @test backend isa CUDABackend
    end
end
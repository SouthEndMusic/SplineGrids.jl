using SafeTestsets

@safetestset "Utilities" include("test_utils.jl")
@safetestset "KnotVector" include("test_knot_vector.jl")
@safetestset "SplineDimension" include("test_spline_dimension.jl")
@safetestset "SplineGrid" include("test_spline_grid.jl")
@safetestset "NURBSGrid" include("test_nurbs_grid.jl")
@safetestset "RefinementMatrix" include("test_refinement_matrix.jl")
@safetestset "Refinement" include("test_refinement.jl")
@safetestset "Local refinement" include("test_local_refinement.jl")
@safetestset "Plotting" include("test_plotting.jl")
@safetestset "LinearMapsExt" include("test_LinearMapsExt.jl")
@safetestset "EnzymeExt" include("test_EnzymeExt.jl")
@safetestset "Aqua" include("aqua.jl")

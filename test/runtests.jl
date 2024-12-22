using SafeTestsets

@safetestset "Utilities" include("test_utils.jl")
@safetestset "KnotVector" include("test_knot_vector.jl")
@safetestset "SplineDimension" include("test_spline_dimension.jl")
@safetestset "SplineGrid" include("test_spline_grid.jl")
@safetestset "LinearMapsExt" include("test_LinearMapsExt.jl")
@safetestset "EnzymeExt" include("test_EnzymeExt.jl")
@safetestset "Aqua" include("aqua.jl")

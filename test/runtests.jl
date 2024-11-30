using SafeTestsets

@safetestset "KnotVector" include("test_knot_vector.jl")
@safetestset "SplineDimension" include("test_spline_dimension.jl")

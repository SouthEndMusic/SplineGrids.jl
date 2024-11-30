using SafeTestsets

@safetestset "SplineDimension" include("test_spline_dimension.jl")
@safetestset "KnotVector" include("test_knot_vector.jl")
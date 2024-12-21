using SafeTestsets

@safetestset "KnotVector" include("test_knot_vector.jl")
@safetestset "SplineDimension" include("test_spline_dimension.jl")
@safetestset "SplineGrid" include("test_spline_grid.jl")
@safetestset "Aqua" include("aqua.jl")

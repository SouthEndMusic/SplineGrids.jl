using SafeTestsets

@safetestset "Utilities" include("test_utils.jl")
@safetestset "KnotVector" include("test_knot_vector.jl")
@safetestset "SplineDimension" include("test_spline_dimension.jl")
@safetestset "SplineGrid" include("test_spline_grid.jl")
@safetestset "NURBSGrid" include("test_nurbs_grid.jl")
@safetestset "Plotting" include("test_plot_rec.jl")
@safetestset "LinearMapsExt" include("test_LinearMapsExt.jl")
@safetestset "EnzymeExt" include("test_EnzymeExt.jl")
@safetestset "Aqua" include("aqua.jl")

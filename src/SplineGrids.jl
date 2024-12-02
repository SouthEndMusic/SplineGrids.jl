module SplineGrids

using KernelAbstractions
using RecipesBase

include("knot_vector.jl")
include("spline_dimension.jl")
include("spline_grid.jl")
include("plot_rec.jl")

export KnotVector, SplineDimension, SplineGrid, decompress, evaluate!

end # module SplineGrids

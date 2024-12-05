module SplineGrids

using KernelAbstractions
using RecipesBase
using PrettyTables
using Subscripts

include("knot_vector.jl")
include("spline_dimension.jl")
include("spline_grid.jl")
include("plot_rec.jl")
include("utils.jl")

export KnotVector, SplineDimension, SplineGrid, decompress, evaluate!

end # module SplineGrids

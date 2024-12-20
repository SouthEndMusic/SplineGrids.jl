module SplineGrids

using Atomix
using KernelAbstractions
using RecipesBase
using PrettyTables
using Subscripts

abstract type AbstractSplineDimension end
abstract type AbstractSplineGrid{Nin, Nout} end

include("utils.jl")
include("knot_vector.jl")
include("spline_dimension.jl")
include("spline_grid.jl")
include("plot_rec.jl")
include("validation.jl")

export KnotVector, SplineDimension, SplineGrid, decompress, evaluate!, evaluate_adjoint!

end # module SplineGrids

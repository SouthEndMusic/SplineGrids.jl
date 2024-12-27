module SplineGrids

using Atomix
using KernelAbstractions

# Needed for pretty printing
using PrettyTables
using Subscripts

# Needed for plotting
using RecipesBase

# Needed for refinement
using SparseArrays
using LinearAlgebra

abstract type AbstractSplineDimension end
abstract type AbstractSplineGrid{Nin, Nout, HasWeights} end
const AbstractNURBSGrid = AbstractSplineGrid{Nin, Nout, true} where {Nin, Nout}

include("utils.jl")
include("knot_vector.jl")
include("spline_dimension.jl")
include("refinement.jl")
include("spline_grid.jl")
include("nurbs_grid.jl")
include("plot_rec.jl")
include("validation.jl")

export KnotVector, SplineDimension, SplineGrid, NURBSGrid, decompress, evaluate!,
       evaluate_adjoint!, insert_knot!, refine!

end # module SplineGrids

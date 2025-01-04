module SplineGrids

# Needed for Kernel computations and memory handling
using Adapt
using Atomix
using KernelAbstractions

# Needed for pretty printing
using PrettyTables
using Subscripts

# Needed for plotting
using RecipesBase

# Needed for updating objects
using Accessors

# Type parameters:
# - Nin: Number of input parameters
# - Nout: number of output parameters
# - HasWeights: Whether the spline has weights and thus is a NURBS
# - T: The float type
abstract type AbstractKnotVector{T <: AbstractFloat} end
abstract type AbstractSplineDimension{T <: AbstractFloat} end
abstract type AbstractSplineGrid{Nin, Nout, HasWeights, T <: AbstractFloat} end
const AbstractNURBSGrid = AbstractSplineGrid{Nin, Nout, true, T} where {Nin, Nout, T}

abstract type AbstractRefinementMatrix{Tv, Ti <: Integer} end

include("util_kernels.jl")
include("utils.jl")
include("knot_vector.jl")
include("spline_dimension.jl")
include("refinement_matrix.jl")
include("refinement.jl")
include("spline_grid.jl")
include("nurbs_grid.jl")
include("plot_rec.jl")
include("validation.jl")

export KnotVector, SplineDimension, SplineGrid, NURBSGrid, decompress, evaluate!,
       evaluate_adjoint!, insert_knot, refine, RefinementMatrix, rmeye, mult!

end # module SplineGrids

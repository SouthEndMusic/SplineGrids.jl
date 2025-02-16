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

# Needed for allocating static arrays in kernels
using StaticArraysCore

# Needed for updating objects
using ConstructionBase

# Needed for handling locally refined control point data
using LazyArrays

# Type parameters:
# - Nin: Number of input parameters
# - Nout: number of output parameters
# - HasWeights: Whether the spline has weights and thus is a NURBS
# - Tv: The floating punt number type
# - Ti: The integer type
abstract type AbstractKnotVector{backend, Tv <: AbstractFloat, Ti <: Integer} end
abstract type AbstractSplineDimension{
    backend, Tv <: AbstractFloat, Ti <: Integer} end
abstract type AbstractSplineGrid{
    Nin, Nout, backend, Tv <: AbstractFloat, Ti <: Integer, HasWeights} end
abstract type AbstractRefinementMatrix{backend, Tv, Ti <: Integer} end
abstract type AbstractControlPoints{Nin, Nout, backend, Tv <: AbstractFloat} end

const AbstractControlPointArray{Nin, Nout, Tv, backend} = Union{
    SplineGrids.AbstractControlPoints{Nin, Nout, backend, Tv},
    AbstractArray{Tv}
} where {Nin, Nout, Tv, backend}

const AbstractNURBSGrid{Nin, Nout, backend, Tv, Ti} = AbstractSplineGrid{
    Nin, Nout, backend, Tv, Ti, true} where {Nin, Nout, backend, Tv, Ti}

include("util_kernels.jl")
include("utils.jl")
include("knot_vector.jl")
include("spline_dimension.jl")
include("refinement_matrix.jl")
include("refinement.jl")
include("control_points.jl")
include("spline_grid.jl")
include("adjoint.jl")
include("plot_rec.jl")
include("validation.jl")

export activate_local_control_point_range!,
       activate_local_refinement!,
       add_default_local_refinement,
       deactivate_overwritten_control_points!,
       decompress,
       DefaultControlPoints,
       error_informed_local_refinement!,
       evaluate_adjoint!,
       evaluate!,
       get_n_control_points,
       insert_knot,
       KnotVector,
       LocallyRefinedControlPoints,
       LocalRefinement,
       mult!,
       NURBSGrid,
       plot_basis,
       plot_basis!,
       refine,
       RefinementMatrix,
       rmeye,
       SplineDimension,
       SplineGrid

# Define names for SplineGridsMakieExt
function plot_basis end
function plot_basis! end

end # module SplineGrids

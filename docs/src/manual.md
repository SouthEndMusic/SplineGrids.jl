# Constructors

```@docs
KnotVector(::Integer, ::Integer)
SplineDimension(::Integer, ::Integer, ::Integer)
SplineGrid(::NTuple{Nin, <:SplineDimension{Tv, Ti}}, ::Integer) where {Nin, Tv, Ti}
NURBSGrid
rmeye
RefinementMatrix(::SplineDimension{Tv, Ti}, ::Integer, ::Any) where {Tv, Ti <: Integer}
```

# Evaluation

```@docs
evaluate!(::SplineDimension)
evaluate!(::SplineGrids.AbstractSplineGrid{Nin, Nout, false}) where {Nin, Nout}
evaluate!(::SplineGrids.AbstractNURBSGrid{Nin}) where {Nin}
evaluate_adjoint!
mult!
evaluate!(::LocallyRefinedControlPoints)
```

# Geometric operations

```@docs
insert_knot(::KnotVector, ::AbstractFloat)
insert_knot(::SplineDimension{Tv, Ti}, ::AbstractFloat) where {Tv, Ti}
insert_knot(::A, ::Integer, ::AbstractFloat) where {A <: SplineGrids.AbstractSplineGrid}
refine(::SplineDimension{Tv, Ti}) where {Tv, Ti}
refine(::A, ::Integer) where {A <: SplineGrids.AbstractSplineGrid}
refine(::SplineGrids.AbstractSplineGrid{Nin, Nout, false, Tv, Ti}, ::SplineDimension{Tv, Ti}, ::Integer, ::RefinementMatrix{Tv, Ti}) where {Nin, Nout, Tv, Ti}
activate_local_refinement!(
        ::LocallyRefinedControlPoints{Nin, Nout, Tv, Ti},
        ::AbstractMatrix{Ti}) where {Nin, Nout, Tv, Ti}
activate_local_refinement!(::SplineGrids.AbstractSplineGrid, args...)
activate_local_control_point_range!
```

# Structs

```@docs
KnotVector
SplineDimension
RefinementMatrix
DefaultControlPoints
LocalRefinement
LocallyRefinedControlPoints
SplineGrid
```

# Utility functions

```@docs
decompress
get_n_control_points
add_default_local_refinement
```

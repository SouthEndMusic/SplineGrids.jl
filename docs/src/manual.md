# Constructors

```@docs
KnotVector(::Integer, ::Integer)
SplineDimension(::Integer, ::Integer, ::Integer)
SplineGrid(::NTuple{Nin, <:SplineDimension{Tv, Ti}}, ::Integer) where {Nin, Tv, Ti}
rmeye(::Integer)
RefinementMatrix(::SplineDimension{Tv, Ti}, ::Integer, ::Any) where {Tv, Ti <: Integer}
```

# Evaluation

```@docs
evaluate!(::SplineDimension)
evaluate!(::SplineGrids.AbstractSplineGrid{Nin, Nout, false}) where {Nin, Nout}
evaluate!(SplineGrids.AbstractNURBSGrid{Nin}) where {Nin}
evaluate_adjoint!(::SplineGrid{Nin}) where {Nin}
mult!(::A, ::RefinementMatrix{Tv}, ::A, ::Integer) where {Tv, A <: AbstractArray{Tv}}
```

# Geometric operations

```@docs
insert_knot(::KnotVector, ::AbstractFloat)
insert_knot(::SplineDimension{Tv, Ti}, ::AbstractFloat) where {Tv, Ti}
insert_knot(::A, ::Integer, ::AbstractFloat) where {A <: SplineGrids.AbstractSplineGrid}
refine(::SplineDimension{Tv, Ti}) where {Tv, Ti}
refine(::A, ::Integer) where {A <: SplineGrids.AbstractSplineGrid}
refine(::SplineGrids.AbstractSplineGrid{Nin, Nout, false, Tv, Ti}, ::SplineDimension{Tv, Ti}, ::Integer, ::RefinementMatrix{Tv, Ti}) where {Nin, Nout, Tv, Ti}
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
refine_all_dimensions_locally
```
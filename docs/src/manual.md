# Constructors

```@docs
KnotVector(::Integer, ::Integer)
SplineDimension(::Integer, ::Integer, ::Integer)
SplineGrid(::NTuple{Nin, <:SplineDimension}, ::Integer) where Nin
rmeye(::Integer)
RefinementMatrix(::SplineGrids.AbstractSplineDimension{Tv}, ::Ti, ::Any) where {Tv, Ti <: Integer}
```

# Evaluation

```@docs
evaluate!(::SplineDimension)
evaluate!(::SplineGrids.AbstractSplineGrid{Nin, Nout}) where {Nin, Nout}
evaluate!(::SplineGrids.AbstractSplineGrid{Nin, Nout, true})  where {Nin, Nout}
evaluate_adjoint!(::SplineGrids.AbstractSplineGrid{Nin, Nout}) where{Nin, Nout}
mult!(::A, ::SplineGrids.AbstractRefinementMatrix{Tv}, ::A, ::Integer) where {Tv, A <: AbstractArray{Tv}}
```

# Geometric operations

```@docs
insert_knot(::KnotVector, ::AbstractFloat)
insert_knot(::SplineDimension, ::AbstractFloat)
insert_knot(::SplineGrids.AbstractSplineGrid{Nin, Nout, HasWeights, T}, ::Integer, ::AbstractFloat) where {Nin, Nout, HasWeights, T}
refine(::SplineGrids.AbstractSplineDimension{Tv}) where {Tv}
refine(::SplineGrids.AbstractSplineGrid{Nin, Nout, HasWeights, T}, ::Integer) where {Nin, Nout, HasWeights, T}
refine(::SplineGrids.AbstractSplineGrid{Nin, Nout, false, T}, ::SplineGrids.AbstractSplineDimension{T}, ::Integer, ::RefinementMatrix) where {Nin, Nout, T}
```

# Structs

```@docs
KnotVector
SplineDimension
SplineGrid
RefinementMatrix
```

# Utility functions

```@docs
decompress
```
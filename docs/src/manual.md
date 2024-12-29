# Constructors

```@docs
KnotVector(::Integer, ::Integer)
SplineDimension(::Integer, ::Integer, ::Integer)
SplineGrid(::NTuple{Nin, <:SplineDimension}, ::Integer) where Nin
```

# Evaluation

```@docs
evaluate!(::SplineDimension)
evaluate!(::SplineGrids.AbstractSplineGrid{Nin, Nout}) where {Nin, Nout}
evaluate!(::SplineGrids.AbstractSplineGrid{Nin, Nout, true})  where {Nin, Nout}
evaluate_adjoint!(::SplineGrids.AbstractSplineGrid{Nin, Nout}) where{Nin, Nout}
```

# Geometric operations

```@docs
insert_knot!(::KnotVector, ::AbstractFloat)
insert_knot!(::SplineDimension, ::AbstractFloat)
insert_knot!(::SplineGrids.AbstractSplineGrid{Nin, Nout}, ::Integer, ::AbstractFloat) where {Nin, Nout}
refine!(::SplineDimension)
refine!(::SplineGrids.AbstractSplineGrid, ::Integer)
```

# Structs

```@docs
KnotVector
SplineDimension
SplineGrid
```

# Utility functions

```@docs
decompress
```
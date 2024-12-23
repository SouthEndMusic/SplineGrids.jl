# Constructors

```@docs
KnotVector(::Integer, ::Integer)
SplineDimension(::Integer, ::Integer, ::Integer)
SplineGrid(::NTuple{N_in, <:SplineDimension}, ::Integer) where N_in
```

# Evaluation

```@docs
evaluate!(::SplineDimension)
evaluate!(::SplineGrids.AbstractSplineGrid{Nin, Nout}) where {Nin, Nout}
evaluate!(SplineGrids.AbstractSplineGrid{Nin, Nout, true})  where {Nin, Nout}
evaluate_adjoint!(::SplineGrids.AbstractSplineGrid{Nin, Nout}) where{Nin, Nout}
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
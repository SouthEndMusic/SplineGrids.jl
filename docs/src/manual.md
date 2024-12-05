# Constructors

```@docs
KnotVector(::Integer, ::Integer)
SplineDimension(::Integer, ::Integer, ::Integer)
SplineGrid(::NTuple{N_in, <:SplineDimension}, ::Integer) where N_in
```

# Evaluation

```@docs
evaluate!(::SplineDimension)
evaluate!(::SplineGrid)
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
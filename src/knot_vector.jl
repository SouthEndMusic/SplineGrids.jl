"""
    KnotVector(knots, multiplicities)

Defines a knot vector.

## Arguments

  - `knot_values`: The values in the knot vector. Must be strictly increasing.
  - `multiplicities`: The multiplicity of each knot in `knots`.
"""
struct KnotVector{
    backend,
    Tv <: AbstractFloat,
    Ti <: Integer,
    K <: AbstractVector{Tv},
    M <: AbstractVector{Ti}
} <: AbstractKnotVector{backend, Tv, Ti}
    knot_values::K
    multiplicities::M
    knots_all::K
    extent::Tuple{Tv, Tv}
    function KnotVector(knot_values, multiplicities)
        @assert length(knot_values)==length(multiplicities) "knot_values and multiplicities must be of the same length."
        @assert knot_values==sort(knot_values) "knot_values must be sorted."
        # allunique cannot be performed efficiently with a CuArray
        knot_values_cpu = adapt(CPU(), knot_values)
        @assert allunique(knot_values_cpu) "knot_values must be unique."

        # Construct explicit knot vector
        knots_all = similar(knot_values, sum(multiplicities))
        backend = get_backend(knot_values)
        expand_knot_vector_kernel(backend)(
            knots_all,
            knot_values,
            multiplicities,
            ndrange = size(knot_values)
        )
        synchronize(backend)
        new{
            backend,
            eltype(knot_values),
            eltype(multiplicities),
            typeof(knot_values),
            typeof(multiplicities)
        }(
            knot_values,
            multiplicities,
            knots_all,
            (first(knot_values_cpu), last(knot_values_cpu))
        )
    end
end

"""
    KnotVector(
        n_basis_functions::Integer, 
        degree::Integer; 
        extent::Tuple{Number, Number} = (0,1), 
        distribution::Symbol = :equispaced)

Construct a clamped knot vector, i.e. the multiplicity of the first and last knot is degree + 1 and the
other multiplicities are 1.

## Arguments

  - `n_basis_functions`: The number of basis functions that will be defined on this knot vector
  - `degree`: The degree of the basis functions that will be defined on this knot vector

## Keyword Arguments

  - `extent`: A tuple (t_min, t_max) defining the extend of the knot vector
  - `distribution`: The distribution of the internal knots. The options are :equispaced or :random
  - `backend`: The KernelAbstractions backend of the arrays in the object. Defaults to `CPU()`.
  - `float_type`: The type of all floating point arrays. Defaults to `Float32`.
  - `int_type`: The type of all integer arrays. Defaults to `Int32`.
"""
function KnotVector(
        n_basis_functions::Integer,
        degree::Integer;
        extent::Tuple{Number, Number} = (0, 1),
        distribution::Symbol = :equispaced,
        backend::Backend = CPU(),
        float_type::Type{Tv} = Float32,
        int_type::Type{Ti} = Int32
)::KnotVector where {Tv <: AbstractFloat, Ti <: Integer}
    @assert n_basis_functions - degree ≥ 1
    n_knot_values = n_basis_functions - degree + 1

    if distribution == :random
        knot_values = cumsum(rand(float_type, n_knot_values))
        knot_values .-= knot_values[1]
        knot_values ./= knot_values[end] - knot_values[1]
        knot_values .*= extent[2] - extent[1]
        knot_values .+= extent[1]

    elseif distribution == :equispaced
        knot_values = float_type.(
            collect(
            range(extent...,
            length = n_knot_values)
        )
        )
    else
        error("Unsupported knot distribution type $distribution.")
    end

    # For now always create a clamped knot vector.
    multiplicities = ones(int_type, n_knot_values)
    multiplicities[1] = degree + 1
    multiplicities[end] = degree + 1
    KnotVector(
        adapt(backend, knot_values),
        adapt(backend, multiplicities)
    )
end

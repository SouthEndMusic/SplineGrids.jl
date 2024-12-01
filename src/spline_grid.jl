struct SplineGrid{S <: SplineDimension, C <: AbstractArray, W <: Union{AbstractArray, Nothing}, E <: AbstractArray, B <: AbstractArray, N_in}
    spline_dimensions::NTuple{N_in, S}
    control_points::C
    weights::W
    eval::E
    basis_function_products::B
    function SplineGrid(spline_dimensions, control_points, weights, eval, basis_function_products)
        # TODO: Add validation of combination of control points, weights, and basis functions
        N_in = length(spline_dimensions)
        new{eltype(spline_dimensions), typeof(control_points), typeof(weights), typeof(eval), typeof(basis_function_products), N_in}(spline_dimensions, control_points, weights, eval, basis_function_products)
    end
end

# TODO: Add indexing functionality such that indexing a spline_grid implies indexing spline_grid.eval

function SplineGrid(spline_dimensions::NTuple{N_in, <:SplineDimension}, dim_out::Integer)::SplineGrid where N_in
    # The size of the point grid on which the spline is evaluated
    size_eval_grid = (length(spline_dimension.sample_points) for spline_dimension in spline_dimensions)
    # The size of the grid of control points
    size_cp_grid = get_n_basis_functions.(spline_dimensions)
    # The actual control points
    control_points = zeros(size_cp_grid..., dim_out)
    # Preallocated memory for basis function product evaluation
    basis_function_products = zeros(size_eval_grid...)
    # Preallocated memory for grid evaluation of the spline
    eval = zeros(size_eval_grid..., dim_out)
    # NURBS are not supported yet
    weights = nothing
    SplineGrid(spline_dimensions, control_points, weights, eval, basis_function_products)
end

cp_kernel_size(spline_grid::SplineGrid) = Tuple(spline_dim.degree + 1 for spline_dim in spline_grid.spline_dimensions)

# With thanks to Michael Abbott
function outer!(A::AbstractArray{T,N}, vs::Vararg{SubArray,N}) where {T,N}
    vecs = ntuple(n -> reshape(vs[n], ntuple(Returns(1), n-1)..., :), N)
    broadcast!(*, A, vecs...)
end

function evaluate!(spline_grid::SplineGrid)::Nothing
    (; basis_function_products, eval, spline_dimensions, control_points) = spline_grid
    eval .= 0

    for I in CartesianIndices(cp_kernel_size(spline_grid))
        outer!(basis_function_products, (view(spline_dim.eval, :, i) for (i, spline_dim) in zip(Tuple(I), spline_dimensions))...)
        for J in CartesianIndices(basis_function_products)
            cp_indices = (spline_dim.sample_indices[j] - spline_dim.degree - 1 + i for (i, j, spline_dim) in zip(Tuple(I), Tuple(J), spline_dimensions))
            control_point = view(control_points, cp_indices..., :)
            @. eval[J, :] += basis_function_products[J] * control_point
        end
    end
    return nothing
end
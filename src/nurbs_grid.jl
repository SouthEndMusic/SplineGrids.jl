function NURBSGrid(spline_dimensions::NTuple{Nin, <:SplineDimension},
        Nout::Integer)::AbstractSplineGrid{Nin} where {Nin}
    # The size of the point grid on which the spline is evaluated
    size_eval_grid = get_sample_grid_size(spline_dimensions)
    # The size of the grid of control points
    size_cp_grid = get_n_basis_functions.(spline_dimensions)
    # The control points
    control_points = zeros(size_cp_grid..., Nout)
    set_unit_cp_grid!(control_points)
    # The weights and denominator
    weights = ones(size_cp_grid...)
    denominator = zeros(size_eval_grid...)
    # Preallocated memory for basis function product evaluation
    basis_function_products = zeros(size_eval_grid...)
    # Preallocated memory for grid evaluation of the spline
    eval = zeros(size_eval_grid..., Nout)
    # Linear indices for control points per global sample point
    sample_indices = get_global_sample_indices(spline_dimensions, control_points)
    SplineGrid(
        spline_dimensions,
        control_points,
        eval,
        sample_indices,
        basis_function_products;
        weights,
        denominator
    )
end

# Support inputting a single spline dimension instead of a tuple
function NURBSGrid(spline_dimension::SplineDimension, args...; kwargs...)
    NURBSGrid((spline_dimension,), args...; kwargs...)
end

@kernel function nurbs_muladd_kernel(
        eval,
        denominator,
        @Const(basis_function_products),
        @Const(weights),
        @Const(control_points),
        @Const(sample_indices),
        offset,
        divide
)
    # Index of the global sample point
    J = @index(Global, Cartesian)

    # Output dimensionality
    Nout = size(control_points)[end]

    # The total number of control points
    cp_grid_size = prod(size(control_points)[1:(end - 1)]) # Could be done outside kernel

    # The product of basis functions for the current sample point
    # and the current control point kernel location
    basis_function_product = basis_function_products[J]

    # The linear index of the required control point
    lin_cp_index_base = sample_indices[J] + offset

    # The weight of the required control_point
    weight = weights[lin_cp_index_base + cp_grid_size]

    # The denominator contribution of the required control point
    denominator_contrib = basis_function_product * weight
    denominator[J] += denominator_contrib

    # The enumerator contribution of the required control point
    for dim_out in 1:Nout
        lin_cp_idx = lin_cp_index_base + dim_out * cp_grid_size
        eval[J, dim_out] += denominator_contrib * control_points[lin_cp_idx]
    end

    # Divide the enumerator by the denominator
    if divide
        denominator_J = denominator[J]
        for dim_out in 1:Nout
            eval[J, dim_out] /= denominator_J
        end
    end
end

"""
    evaluate!(spline_grid::AbstractNURBSGrid;
        control_points::AbstractArray = spline_grid.control_points,
        weights::AbstractArray = spline_grid.weights,
        eval::AbstractArray = spline_grid.eval)

Evaluate the NURBS grid, that is: take the evaluated basis functions for each sample point
for each SplineDimension, and compute the output grid on each sample point combination
as a weighted linear combination of control points with rational functions as coefficients.

Uses the `control_points`, `weights` and `eval` arrays from the `spline_grid` by default,
but different arrays can be specified as a convenience for optimization algorithms.

NOTE: At the moment computing derivatives of NURBS grids is not supported.
"""
function evaluate!(spline_grid::AbstractNURBSGrid;
        control_points::AbstractArray = spline_grid.control_points,
        weights::AbstractArray = spline_grid.weights,
        eval::AbstractArray = spline_grid.eval
)::Nothing
    (; basis_function_products, spline_dimensions, sample_indices, denominator) = spline_grid
    eval .= 0
    denominator .= 0

    @assert size(control_points) == size(spline_grid.control_points)
    @assert size(eval) == size(spline_grid.eval)

    control_point_kernel_size = get_cp_kernel_size(spline_dimensions)
    backend = get_backend(eval)
    kernel! = nurbs_muladd_kernel(backend)

    # Loop over the positions in the kernel of control points each spline evaluation depends on
    for I in CartesianIndices(control_point_kernel_size)
        # Compute basis function products as an outer product of the basis function values per 
        # spline dimension
        outer!(
            basis_function_products,
            (
                view(spline_dim.eval, :, i, 1) for
            (i, spline_dim) in zip(Tuple(I), spline_dimensions))...
        )

        # Divide the enumerator by the denominator after processing
        # the last control point in the control point kernel
        divide = (Tuple(I) == control_point_kernel_size)
        offset = get_offset(size(control_points), Tuple(I))
        kernel!(
            eval,
            denominator,
            basis_function_products,
            weights,
            control_points,
            sample_indices,
            offset,
            divide,
            ndrange = size(sample_indices))
        synchronize(backend)
    end
    return nothing
end
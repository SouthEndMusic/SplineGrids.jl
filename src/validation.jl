function validate_spline_grid(
        spline_dimensions,
        control_points,
        denominator,
        weights,
        eval,
        sample_indices,
        basis_function_products
)::Nothing
    errors = false

    # Number of input dimensions
    Nin_sd = length(spline_dimensions)
    Nin_cp = ndims(control_points) - 1
    Nin_eval = ndims(eval) - 1
    Nin_si = ndims(sample_indices)
    Nin_bfp = ndims(basis_function_products)

    if isnothing(weights)
        Nins = (Nin_sd, Nin_cp, Nin_eval, Nin_si, Nin_bfp)
        if !allequal(Nins)
            @error "The number of input dimensions from the spline dimensions, control_points, eval, sample indices and basis function products must agree, got $Nins respectively."
            errors = true
        end
    else
        Nin_weights = ndims(weights)
        Nin_denom = ndims(weights)
        Nins = (Nin_sd, Nin_cp, Nin_eval, Nin_si, Nin_bfp, Nin_weights, Nin_denom)
        if !allequal(Nins)
            @error "The number of input dimensions from the spline dimensions, control_points, eval, sample indices, basis function products, weights and denominator must agree, got $Nins respectively."
            errors = true
        end
    end

    # Number of output dimensions
    Nout_cp = size(control_points)[end]
    Nout_eval = size(eval)[end]

    if Nout_cp != Nout_eval
        @error "The number of output dimensions from the control points $Nout_cp and eval $Nout_eval must agree."
        errors = true
    end

    # Control point grid size
    cp_grid_size_sd = get_control_point_grid_size(spline_dimensions)
    cp_grid_size_cp = size(control_points)[1:(end - 1)]

    if cp_grid_size_sd != cp_grid_size_cp
        @error "The control point grid sizes from the spline dimensions $cp_grid_size_sd and the control points $cp_grid_size_cp must agree."
        errors = true
    end

    # Evaluation grid size
    eval_grid_size_si = size(sample_indices)
    eval_grid_size_bfp = size(basis_function_products)
    eval_grid_size_eval = size(eval)[1:(end - 1)]

    if isnothing(weights)
        eval_grid_sizes = (
            eval_grid_size_si,
            eval_grid_size_bfp,
            eval_grid_size_eval
        )
        if !allequal(eval_grid_sizes)
            @error "The evaluation grid sizes from the sample indices, basis function products and eval must agree, got $eval_grid_sizes respectively."
            errors = true
        end
    else
        eval_grid_size_denom = size(denominator)
        eval_grid_sizes = (
            eval_grid_size_si,
            eval_grid_size_bfp,
            eval_grid_size_eval,
            eval_grid_size_denom
        )
        if !allequal(eval_grid_sizes)
            @error "The evaluation grid sizes from the sample indices, basis function products, eval and denominator must agree, got $eval_grid_sizes respectively."
            errors = true
        end
    end

    if errors
        error("Errors encountered when validating spline grid constructor inputs.")
    end
end

function validate_partial_derivatives(
        spline_dimensions::NTuple{Nin, <:SplineDimension},
        derivative_order::NTuple{
            Nin, <:Integer}
)::Nothing where {Nin}
    errors = false
    for (i, (spline_dimension, derivative_order_)) in enumerate(zip(
        spline_dimensions, derivative_order))
        (; max_derivative_order) = spline_dimension

        if !(0 ≤ derivative_order_ ≤ max_derivative_order)
            errors = true
            @error "The maximum derivative order available for spline dimension $i is $max_derivative_order, got $derivative_order_."
        end
    end

    if errors
        error("Invalid derivative order(s) supplied. If you want to evaluate (higher order) derivatives, specify this at construction as SplineDimension(...; max_derivative_order).")
    end
    return nothing
end

function validate_partial_derivatives(
        spline_grid::AbstractSplineGrid{Nin, Nout, false},
        derivative_order::NTuple{Nin, <:Integer}
)::Nothing where {Nin, Nout}
    validate_partial_derivatives(spline_grid.spline_dimensions, derivative_order)
end

function validate_partial_derivatives(
        ::AbstractNURBSGrid,
        derivative_order::NTuple{Nin, <:Integer}
)::Nothing where {Nin}
    if any(!iszero, derivative_order)
        error("Computing derivatives of NURBS is currently not supported.")
    end
end
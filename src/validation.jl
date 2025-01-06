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
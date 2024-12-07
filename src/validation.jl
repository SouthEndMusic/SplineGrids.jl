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
            @error "The maximum derivative order available for spline dimension $i is $max_derivative_order, got $derivative_order_. If you want to evaluate higher order derivatives, specify this at construction as SplineDimension(...; max_derivative_order)."
        end
    end

    if errors
        error("Invalid derivative order(s) supplied.")
    end
    return nothing
end
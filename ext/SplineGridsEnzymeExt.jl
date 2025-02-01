module SplineGridsEnzymeExt
using SplineGrids
using Enzyme
using .EnzymeRules

function Enzyme.make_zero!(spline_grid::SplineGrid)::Nothing
    make_zero!(spline_grid.eval)
    for spline_dimension in spline_grid.spline_dimensions
        make_zero!(spline_dimension.eval)
    end
    return nothing
end

function EnzymeRules.augmented_primal(
        config::RevConfigWidth{1},
        ::Const{typeof(evaluate!)},
        ::Type{RT},
        spline_grid::Duplicated{<:SplineGrid},
        control_points::Duplicated{<:SplineGrids.AbstractControlPointArray};
        kwargs...
) where {RT}
    println("yeet")
    evaluate!(spline_grid.val, control_points.val; kwargs...)
    primal = if needs_primal(config)
        spline_grid.val
    else
        nothing
    end
    shadow = if needs_shadow(config)
        spline_grid.dval
    else
        nothing
    end
    EnzymeRules.AugmentedReturn(primal, shadow, kwargs)
end

function EnzymeRules.reverse(
        ::RevConfigWidth{1},
        ::Const{typeof(evaluate!)},
        ::Type{RT},
        kwargs,
        spline_grid::Duplicated{<:SplineGrid},
        control_points::Duplicated{<:SplineGrids.AbstractControlPointArray}
) where {RT}
    println("yaat")
    evaluate_adjoint!(spline_grid.dval; control_points, kwargs...)
    (nothing, nothing)
end

end # Module SplineGridsEnzymeExt

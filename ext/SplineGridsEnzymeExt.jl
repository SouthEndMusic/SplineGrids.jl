module SplineGridsEnzymeExt
using SplineGrids
using Enzyme
using .EnzymeRules
using ConstructionBase

function Enzyme.make_zero(spline_grid::SplineGrid)
    setproperties(
        spline_grid;
        eval = make_zero(spline_grid.eval),
        control_points = make_zero(spline_grid.control_points)
    )
end

function Enzyme.make_zero!(spline_grid::SplineGrid)::Nothing
    make_zero!(spline_grid.eval)
    make_zero!(spline_grid.control_points)
    return nothing
end

function EnzymeRules.augmented_primal(
        config::RevConfigWidth{1},
        ::Const{typeof(evaluate!)},
        ::Type{Const{Nothing}},
        spline_grid::MixedDuplicated{<:SplineGrid},
        control_points::Duplicated{<:SplineGrids.AbstractControlPointArray};
        kwargs...
)
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
        ::Type{Const{Nothing}},
        kwargs,
        spline_grid::MixedDuplicated{<:SplineGrid},
        control_points::Duplicated{<:SplineGrids.AbstractControlPointArray}
)
    evaluate_adjoint!(spline_grid.dval[]; control_points = control_points.dval, kwargs...)
    make_zero!(spline_grid.dval[])
    (nothing, nothing)
end

end # Module SplineGridsEnzymeExt

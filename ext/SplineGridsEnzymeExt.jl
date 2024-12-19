module SplineGridsEnzymeExt
using SplineGrids
using Enzyme

function Enzyme.make_zero!(spline_grid::SplineGrid)::Nothing
    make_zero!(spline_grid.eval)
    make_zero!(spline_grid.basis_function_products)
    for spline_dimension in spline_grid.spline_dimensions
        make_zero!(spline_dimension.eval)
    end
    return nothing
end

end # Module SplineGridsEnzymeExt
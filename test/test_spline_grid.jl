using SplineGrids

@testset "Spline curve in 1D" begin
    n_sample_points = 100
    for n_basis_functions in 2:10
        for degree in 1:(n_basis_functions - 1)
            spline_dimension = SplineDimension(n_basis_functions, degree, n_sample_points)
            spline_grid = SplineGrid(spline_dimension, 1)
            spline_grid.control_points .= 1
            evaluate!(spline_grid)
            @test all(spline_grid.eval .≈ 1)
        end
    end
end
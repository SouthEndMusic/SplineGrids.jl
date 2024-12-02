using SplineGrids

@testset "Partition of unity" begin
    n_sample_points = 500
    n_basis_functions = 25
    for degree in 0:5
        for distribution in [:equispaced, :random]
            s = SplineDimension(n_basis_functions, degree, n_sample_points; distribution)
            @test all(s.eval .>= 0)
            @test all(sum(s.eval, dims = 2) .≈ 1)
            @test size(s.eval) == (n_sample_points, degree + 1)
        end
    end
end

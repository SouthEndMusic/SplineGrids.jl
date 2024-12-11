using SplineGrids

@testset "Partition of unity" begin
    n_sample_points = 500
    n_basis_functions = 25
    for degree in 0:5
        for distribution in [:equispaced, :random]
            s = SplineDimension(n_basis_functions, degree, n_sample_points; distribution)
            B = s.eval[:, :, 1]
            @test all(B .>= 0)
            @test all(sum(B, dims = 2) .≈ 1)
            @test size(B) == (n_sample_points, degree + 1)
        end
    end
end

@testset "Derivatives" begin
    n_control_points = 10
    degree = 3
    n_sample_points = 5000

    spline_dimension = SplineDimension(
        n_control_points, degree, n_sample_points; max_derivative_order = 2)

    (; sample_points) = spline_dimension
    Δt = diff(sample_points)
    data = decompress(spline_dimension)
    data_deriv = decompress(spline_dimension; derivative_order = 1)
    data_deriv2 = decompress(spline_dimension; derivative_order = 2)

    deriv_fd = diff(data; dims = 1)

    for i in 1:SplineGrids.get_n_basis_functions(spline_dimension)
        deriv_fd[:, i] ./= Δt
    end

    deriv2_fd = diff(deriv_fd; dims = 1)

    for i in 1:SplineGrids.get_n_basis_functions(spline_dimension)
        deriv2_fd[:, i] ./= Δt[2:end]
    end

    using Test
    @test data_deriv[1:(end - 1), :]≈deriv_fd rtol=1e-2
    @test data_deriv2[3:end, :]≈deriv2_fd rtol=1e-2
end

using SplineGrids

@testset "Validation" begin
    knot_values = collect(1:5)
    multiplicities = collect(1:3)
    @test_throws "knot_values and multiplicities must be of the same length." KnotVector(knot_values, multiplicities)

    knot_values = [3,2,1]
    multiplicities = collect(1:3)
    @test_throws "knot_values must be sorted." KnotVector(knot_values, multiplicities)

    knot_values = [1,1,2]
    multiplicities = collect(1:3)
    @test_throws "knot_values must be unique." KnotVector(knot_values, multiplicities)
end

@testset "Construction" begin
    knot_vector_1 = KnotVector(5, 2)
    @test knot_vector_1.knot_values ≈ 0.0:(1/3):1.0
    @test knot_vector_1.multiplicities == [3, 1, 1, 3]
    @test knot_vector_1.knots_all ≈ [0.0, 0.0, 0.0, 1/3, 2/3, 1.0, 1.0, 1.0]

    knot_vector_2 = KnotVector(5, 2, extent = (5.0, 7.0))
    @test knot_vector_2.knot_values ≈ @. 2 * knot_vector_1.knot_values + 5.0 
    @test knot_vector_2.multiplicities == knot_vector_1.multiplicities 
    @test knot_vector_2.knots_all ≈ @. 2 * knot_vector_1.knots_all + 5.0 

    # Just to see whether it passes validation
    knot_vector_3 = KnotVector(10, 3; extent = (4, 8), distribution = :random)

    @test_throws AssertionError KnotVector(5, 10)
end
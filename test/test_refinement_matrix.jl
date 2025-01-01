using SplineGrids
using KernelAbstractions
using LinearAlgebra
using BandedMatrices
using Adapt
using Random

Random.seed!(1)

if "--gpu_backend" ∈ ARGS
    backend = CUDABackend()
else
    backend = CPU()
end

@testset "Unit matrix" begin
    N = 100
    M = RefinementMatrix(collect(Float64.(I(N))))
    M = adapt(backend, M)
    M_squared = M * M
    @test M_squared isa RefinementMatrix
    @test M_squared.m == M.m
    @test M_squared.n == M.n
    @test M_squared.row_pointer == M.row_pointer
    @test M_squared.column_start == M.column_start
    @test M_squared.nzval == M.nzval
    @test M_squared == M
    @test size(M_squared) == (N, N)
    @test length(M_squared) == N^2
end

@testset "Diagonal matrices" begin
    A = RefinementMatrix(diagm(rand(100)))
    A = adapt(backend, A)

    B = RefinementMatrix(diagm(rand(100)))
    B = adapt(backend, B)

    C = A * B
    @test C.nzval == A.nzval .* B.nzval
end

@testset "Banded matrices" begin
    Random.seed!(1)

    A = RefinementMatrix(collect(brand(100, 100, 4, 3)))
    A = adapt(backend, A)

    B = RefinementMatrix(collect(brand(100, 200, 3, 5)))
    B = adapt(backend, B)

    C = A * B
    @test collect(C) ≈ collect(A) * collect(B)
end

@testset "Empty columns" begin
    M = zeros(2, 50)
    M[1, 1:20] .= 1
    M[2, 39:50] .= 1
    @test_throws "Invalid rows: [2]." RefinementMatrix(M)
end

@testset "Empty rows" begin
    M = zeros(3, 50)
    M[1, 1:25] .= 1
    M[3, 26:50] .= 1
    @test_throws "Invalid rows: [2, 3]." RefinementMatrix(M)
end
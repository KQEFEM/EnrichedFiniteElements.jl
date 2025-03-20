# test/OperatorsTests.jl

using Revise
using Test
using HCubature
using EnrichedFiniteElements

@testset "Operators Tests" begin
    @testset "Reference Square" begin
    # Test Case 1: A = B = C = 0 (Reference Square)
    A = 0.0
    B = 0.0
    C = 0.0
    omega = [0.0, 0.0]
    lower_bounds = [0.0, 0.0, 0.0]
    upper_bounds = [1.0, 1.0, 1.0]

    expected_result1 = zeros(ComplexF64, 3, 3)  # 3x3 matrix of complex zeros
    result1, _ =
        EnrichedFiniteElements.Operators.pDtq(upper_bounds, lower_bounds, A, B, C, omega)

    @test isapprox(real(result1), real(expected_result1), atol = 1e-6)
    @test isapprox(imag(result1), imag(expected_result1), atol = 1e-6)
    end 
    @testset "Complex Valued" begin
    # Test Case 2: Complex A, B, C (Example Complex Values)
    A = 1.0
    B = 0.2
    C = 0.8
    omega = [3.0, -3]
    lower_bounds = [0.0, 0.0, 0.0]
    upper_bounds = [1.0, 1.0, 1.0]
    expected_result2 = ComplexF64[
        -3.7996281927039286e-6 + 0.00018323667891961633im -5.370289355812572e-5 + 0.0009135553137364832im -9.433183249596887e-6 + 0.0003664505797012257im; -5.3702893558125525e-5 + 0.0009135553137364834im -0.0009903075835568681 + 0.007978629788904497im -0.00018308044225269775 + 0.002739782617397888im; -9.433183249596893e-6 + 0.00036645057970122554im -0.00018308044225269813 + 0.002739782617397888im -4.4936269655183895e-5 + 0.0012823678089615273im    ]
    result2, _ =
        EnrichedFiniteElements.Operators.pDtq(upper_bounds, lower_bounds, A, B, C, omega)
    @test isapprox(real(result2), real(expected_result2), atol = 1e-6)
    @test isapprox(imag(result2), imag(expected_result2), atol = 1e-6)

    # Add more test cases as needed (different bounds, omega values, etc.)
end
end 
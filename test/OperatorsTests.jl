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
        t0 = 0
        triangle_area = 1.0

        expected_result1 = zeros(ComplexF64, 3, 3)  # 3x3 matrix of complex zeros
        result1, _ = EnrichedFiniteElements.Operators.pDtq(
            upper_bounds,
            lower_bounds,
            A,
            B,
            C,
            omega,
            t0,
            triangle_area,
        )
        @test isapprox(real(result1), real(expected_result1), atol = 1e-6)
        @test isapprox(imag(result1), imag(expected_result1), atol = 1e-6)
    end
    @testset "Complex Valued" begin
        # Test Case 2: Complex A, B, C (Example Complex Values)
        A = 1.0
        B = 0.2
        C = 0.8
        omega = [3.0, -3]
        t0 = 0.0
        triangle_area = 1.0
        lower_bounds = [0.0, 0.0, 0.0]
        upper_bounds = [1.0, 1.0, 1.0]
        expected_result2 = ComplexF64[
            0.00017700008952752238-4.7550877837710987e-5im 0.0008921740882195902-0.00020369759041829422im 0.0003544907357522281-9.333450902768957e-5im
            0.0008921740882195905-0.00020369759041829466im 0.0079375505383193-0.0012784889010357171im 0.0026818133740891506-0.000589749324279142im
            0.0003544907357522281-9.333450902768963e-5im 0.002681813374089151-0.0005897493242791426im 0.0012438473568946752-0.00031516696929945756im
        ]
        result2, _ = EnrichedFiniteElements.Operators.pDtq(
            upper_bounds,
            lower_bounds,
            A,
            B,
            C,
            omega,
            t0,
            triangle_area,
        )
        println(result2)
        @test isapprox(real(result2), real(expected_result2), atol = 1e-6)
        @test isapprox(imag(result2), imag(expected_result2), atol = 1e-6)

        # Add more test cases as needed (different bounds, omega values, etc.)
    end
end

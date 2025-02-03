# test/OperatorsTests.jl

using Revise
using Test
using HCubature
using EnrichedFiniteElements

@testset "Operators Tests" begin
    # Test Case 1: A = B = C = 0 (Reference Square)
    A = [0.0,0.0]
    B = [0.0,0.0]
    C = [0.0,0.0]
    omega = [0.0,0.0]
    lower_bounds = [0.0, 0.0, 0.0]
    upper_bounds = [1.0, 1.0, 1.0]

    expected_result1 = zeros(ComplexF64, 3, 3)  # 3x3 matrix of complex zeros
        result1,_ = EnrichedFiniteElements.Operators.pDtq(upper_bounds, lower_bounds, A, B, C, omega)

    @test isapprox(real(result1), real(expected_result1), atol=1e-6)
    @test isapprox(imag(result1), imag(expected_result1), atol=1e-6)

    # Test Case 2: Complex A, B, C (Example Complex Values)
    A = [1.0,2] 
    B = [0.2,1] 
    C = [0.8,5] 
    omega = [3.0, - 3]
    expected_result2 = ComplexF64[
        0.010015175735776498 + 0.017124413608400352im -0.013623466942252905 - 0.0034438490238047703im -0.0136926524553957 - 0.003562145634271098im;  
        -0.013623466942252926 - 0.003443849023804772im 0.0316170538330674 + 0.031944221590521515im 0.02489690378178615 + 0.02271969603053285im;  
        -0.0136926524553957 - 0.003562145634271091im 0.02489690378178616 + 0.022719696030532845im 0.029711336001351174 + 0.033217536090198244im   
    ]
    result2,_ = EnrichedFiniteElements.Operators.pDtq(upper_bounds, lower_bounds, A, B, C, omega)

    println(result2)
    @test isapprox(real(result2), real(expected_result2), atol=1e-6)
    @test isapprox(imag(result2), imag(expected_result2), atol=1e-6)

    # Add more test cases as needed (different bounds, omega values, etc.)
end
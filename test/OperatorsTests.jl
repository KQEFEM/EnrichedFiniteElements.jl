# test/OperatorsTests.jl

using Revise
using Test
using HCubature
using EnrichedFiniteElements

@testset "Operators Tests" begin
    # Test Case 1: A = B = C = 0 (Reference Square)
    A = 0.0
    B = 0.0
    C = 0.0
    omega = [0.0,0.0]
    lower_bounds = [0.0, 0.0, 0.0]
    upper_bounds = [1.0, 1.0, 1.0]

    expected_result1 = zeros(ComplexF64, 3, 3)  # 3x3 matrix of complex zeros
        result1,_ = EnrichedFiniteElements.Operators.pDtq(upper_bounds, lower_bounds, A, B, C, omega)

    @test isapprox(real(result1), real(expected_result1), atol=1e-6)
    @test isapprox(imag(result1), imag(expected_result1), atol=1e-6)

    # Test Case 2: Complex A, B, C (Example Complex Values)
    A = 1.0 
    B = 0.2 
    C = 0.8
    omega = [3.0, - 3]
    expected_result2 = ComplexF64[0.0006346012307852497 + 0.021724007406845564im  0.006536407034320188 - 0.010877443789792914im  0.006525621227956039 - 0.011246669259216918im;
    0.006536407034320185 - 0.01087744378979291im  -0.01014101946163837 + 0.0449600976472488im  -0.005867766662986309 + 0.03376669718979576im;
    0.006525621227956048 - 0.011246669259216928im  -0.0058677666629863 + 0.033766697189795775im  -0.0009376097483323095 + 0.04506078838750317im
    ]
    result2,_ = EnrichedFiniteElements.Operators.pDtq(upper_bounds, lower_bounds, A, B, C, omega)

    println(result2)
    @test isapprox(real(result2), real(expected_result2), atol=1e-6)
    @test isapprox(imag(result2), imag(expected_result2), atol=1e-6)

    # Add more test cases as needed (different bounds, omega values, etc.)
end
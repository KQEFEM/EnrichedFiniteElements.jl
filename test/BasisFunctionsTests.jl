module BasisFunctionsTests

using Test
include("../src/BasisFunctions.jl")
using .BasisFunctions  # Assuming BasisFunctions.jl is in the same directory or LOAD_PATH

@testset "BasisFunctions" begin
    @testset "phi" begin
        x = 0.2
        y = 0.3
        phi_val = BasisFunctions.phi(x, y)
        @test size(phi_val) == (3,)  # Check size
        @test phi_val ≈ [0.5; 0.2; 0.3] # Check values (using ≈ for floating-point comparison)
    end

    @testset "p_matrix" begin
        x = 0.2
        y = 0.3
        p_mat = BasisFunctions.p_matrix(x, y)
        @test size(p_mat) == (3, 3)
        phi_val = BasisFunctions.phi(x, y)
        println(phi_val)
        @test p_mat ≈ phi_val * phi_val' # Check against the expected calculation
    end

    @testset "grads_matrix" begin
        grad_matrix_input = [1.0 2.0 3.0; 4.0 5.0 6.0; 7.0 8.0 9.0]
        grad_matrix = BasisFunctions.grads_matrix(grad_matrix_input)
        @test grad_matrix == grad_matrix_input # Check for equality

        # Test for error on incorrect size:
        grad_matrix_bad_size = [1.0 2.0; 3.0 4.0]
        @test_throws ErrorException BasisFunctions.grads_matrix(grad_matrix_bad_size)
    end

    @testset "e_space" begin
        x = 0.2
        y = 0.3
        A = 1.0
        B = 2.0
        C = 0.0
        e_space_val = BasisFunctions.e_space(x, y, A, B, C)
        @test isa(e_space_val, Complex) # Check type
        @test e_space_val ≈ exp(1im * (A * x + B * y + C)) # Check value
    end


    @testset "e_time functions" begin
        t = 0.1
        w = 1.0
        dt = 0.01
        t0 = 0.0
        ww = 2.0

        e_mass = BasisFunctions.e_time_mass(t, w, dt, t0, ww)
        @test isa(e_mass, Complex)

        e_ansatz = BasisFunctions.e_time_ansatz(t, w)
        @test isa(e_ansatz, Complex)

        e_test = BasisFunctions.e_time_test(t, ww)
        @test isa(e_test, Complex)

    end
end
end 
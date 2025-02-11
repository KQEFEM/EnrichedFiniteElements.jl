using Revise
using Test
using LinearAlgebra

using EnrichedFiniteElements

const wavenumber_func = EnrichedFiniteElements.EnrichmentCreator

@testset "wavenumber_creation" begin
    wavenumbers_ansatz = [1 1; 1 -1] # Use Ints from the start
    wavenumbers_test = wavenumbers_ansatz

    ansatz_result, test_result = wavenumber_func.wavenumber_creation(wavenumbers_ansatz, wavenumbers_test)

    expected_result = [
        1.0  1.0 -sqrt(2.0)
        1.0  1.0  sqrt(2.0)
        1.0 -1.0 -sqrt(2.0)
        1.0 -1.0  sqrt(2.0)
    ]

    @test ansatz_result == expected_result
    @test test_result == expected_result

end



@testset "create_wavenumbers" begin
    # Test case 1: Basic test with integer enrichments
    x_enrichments = 1
    y_enrichments = 1

    # Correct way to create the expected array of tuples:
    expected_wavenumbers = [(-1, -1), (0, -1), (1, -1), (-1, 0), (0, 0), (1, 0), (-1, 1), (0, 1), (1, 1)]
    expected_wavenumbers = reshape(expected_wavenumbers, :, 1) # Reshape directly (no collect needed)

    actual_wavenumbers = wavenumber_func.create_wavenumbers(x_enrichments, y_enrichments)
    println(actual_wavenumbers)
    @test actual_wavenumbers == expected_wavenumbers

    # Test case 2: Test with zero enrichments
    x_enrichments = 0
    y_enrichments = 0
    expected_wavenumbers = [(0, 0)]
    expected_wavenumbers = reshape(expected_wavenumbers, :, 1)  # Reshape directly
    actual_wavenumbers = wavenumber_func.create_wavenumbers(x_enrichments, y_enrichments)
    @test actual_wavenumbers == expected_wavenumbers
end

# Unit Tests
@testset "combine_wavenumber_with_all_nodes Tests" begin
    # Test 1: Basic functionality
    mat1 = [1 2; 3 4]  # 2x2 matrix
    mat2 = [5 6; 7 8]  # 2x2 matrix

    expected = [
        ([1, 2], [5, 6]), 
        ([1, 2], [7, 8]), 
        ([3, 4], [5, 6]), 
        ([3, 4], [7, 8])
    ]
    
    result = combine_wavenumber_with_all_nodes(mat1, mat2)
    @test result == expected

    # Test 2: Empty matrix case
    empty_mat = Matrix{Int}(undef, 0, 2)  # Empty 0x2 matrix
    @test combine_wavenumber_with_all_nodes(empty_mat, mat2) == []
    @test combine_wavenumber_with_all_nodes(mat1, empty_mat) == []

    # Test 3: Single-row matrices
    single_row_mat1 = [9 10]
    single_row_mat2 = [11 12]

    @test combine_wavenumber_with_all_nodes(single_row_mat1, single_row_mat2) == [([9, 10], [11, 12])]

    # Test 4: Type stability check
    @test eltype(combine_wavenumber_with_all_nodes(mat1, mat2)) == Tuple{Vector{Int}, Vector{Int}}
end
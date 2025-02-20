using Revise
using LinearAlgebra
using Test

using EnrichedFiniteElements

const wavenumber_func = EnrichedFiniteElements.EnrichmentCreator

@testset "wavenumber_creation" begin
    wavenumbers_ansatz = [1 1; 1 -1] # Use Ints from the start
    wavenumbers_test = wavenumbers_ansatz

    ansatz_result, test_result =
        wavenumber_func.wavenumber_creation(wavenumbers_ansatz, wavenumbers_test)

    expected_result = [
        1.0 1.0 -sqrt(2.0)
        1.0 1.0 sqrt(2.0)
        1.0 -1.0 -sqrt(2.0)
        1.0 -1.0 sqrt(2.0)
    ]

    @test ansatz_result == expected_result
    @test test_result == expected_result

end



@testset "create_wavenumbers" begin
    # Test case 1: Basic test with integer enrichments
    x_enrichments = 1
    y_enrichments = 1

    # Correct way to create the expected array of tuples:
    expected_wavenumbers = [
        -1 0
        -1 1
        0 0
        -1 -1
        1 1
        -1 0
        -1 1
        0 1
        0 1
    ]
    # expected_wavenumbers = reshape(expected_wavenumbers, :, 1) # Reshape directly (no collect needed)
    actual_wavenumbers = wavenumber_func.create_wavenumbers(x_enrichments, y_enrichments)

    @test actual_wavenumbers == expected_wavenumbers

    # Test case 2: Test with zero enrichments
    x_enrichments = 0
    y_enrichments = 0
    expected_wavenumbers = [0 0]
    # expected_wavenumbers = reshape(expected_wavenumbers, :, 1)  # Reshape directly
    actual_wavenumbers = wavenumber_func.create_wavenumbers(x_enrichments, y_enrichments)
    @test actual_wavenumbers == expected_wavenumbers
end

# Unit Tests
@testset "combine_wavenumber_with_all_nodes Tests" begin
    # Test 1: Basic functionality
    mat1 = [1 2; 3 4]  # 2x2 matrix
    mat2 = [5 6; 7 8]  # 2x2 matrix

    expected = reshape(
        [([1, 2], [5, 6]) ([1, 2], [7, 8]); ([3, 4], [5, 6]) ([3, 4], [7, 8])],
        :,
        1,
    )

    result = wavenumber_func.combine_wavenumber_with_all_nodes(mat1, mat2)
    @test result == expected


    # Test 3: Single-row matrices
    single_row_mat1 = [9 10]
    single_row_mat2 = [11 12]
    @test wavenumber_func.combine_wavenumber_with_all_nodes(
        single_row_mat1,
        single_row_mat2,
    ) == [([9, 10], [11, 12]);;]

end

@testset "generate_pairs Function Tests" begin

    vector1_int = [10, 20]
    vector2_int = [1, 2, 3]
    expected_pairs_int = [(10, 1), (20, 1), (10, 2), (20, 2), (10, 3), (20, 3)]
    expected_pairs_int_reshaped = reshape(expected_pairs_int, :, 1)
    actual_pairs_int = wavenumber_func.generate_pairs(vector1_int, vector2_int)
    @test actual_pairs_int == expected_pairs_int_reshaped
end

@testset "Edge Cases - Empty Vectors" begin
    vector1_empty = []
    vector2_non_empty = [1, 2]

    @test_throws MethodError wavenumber_func.generate_pairs(
        vector1_empty,
        vector2_non_empty,
    )  # Test for ArgumentError when vector1 is empty
    @test_throws MethodError wavenumber_func.generate_pairs(
        vector2_non_empty,
        vector1_empty,
    )  # Test for ArgumentError when vector2 is empty
    @test_throws ArgumentError wavenumber_func.generate_pairs(vector1_empty, vector1_empty)      # Test for ArgumentError when both are empty
end

@testset "Connectivity Creation" begin
    # Correct the initialization of 'connectivity' with random integers between 1 and 10
    connectivity = [1 2 3; 4 5 6]  # 2x3 matrix with sequential values

    # Define the vectors
    vector = [1, 2, 3]
    vector_2 = [60, 70]

    # Directly define 'test_data' as a 6×1 matrix of tuples
    test_data = [(1, 60); (2, 60); (3, 60); (1, 70); (2, 70); (3, 70);;]

    # Create the expected Matrix of Tuples by constructing it directly
    test_data_matrix = Matrix{Tuple{Int64,Int64}}(test_data)

    all_pairs = wavenumber_func.generate_pairs(vector, vector_2)

    @test all_pairs == test_data
    result = wavenumber_func.combine_wavenumber_with_all_nodes(all_pairs, connectivity)

    test_data = [
        ([(1, 60)], [1, 2, 3]); ([(2, 60)], [1, 2, 3]); ([(3, 60)], [1, 2, 3]); ([(1, 70)], [1, 2, 3]); ([(2, 70)], [1, 2, 3]); ([(3, 70)], [1, 2, 3]); ([(1, 60)], [4, 5, 6]); ([(2, 60)], [4, 5, 6]); ([(3, 60)], [4, 5, 6]); ([(1, 70)], [4, 5, 6]); ([(2, 70)], [4, 5, 6]); ([(3, 70)], [4, 5, 6]);;
    ]

    @test result == test_data


end

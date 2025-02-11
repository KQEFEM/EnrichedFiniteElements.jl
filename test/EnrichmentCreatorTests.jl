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

@testset "combine_wavenumber_with_all_nodes" begin
    # Smaller example wavenumbers (2 wavenumbers)
    wavenumbers = [(-1, -1), (0, -1)]
    wavenumbers = reshape(wavenumbers, :, 1)
    println(wavenumbers)

    # Smaller example nodes (2 nodes)
    nodes = [[0.0, 0.0], [1.0, 0.0]]
    nodes = reshape(nodes, :, 1)

    # Call the function
    combined_result = wavenumber_func.combine_wavenumber_with_all_nodes(wavenumbers, nodes)
    println(combined_result)
    # Expected result (manually constructed for clarity)
    expected_result = [
        ((-1, -1), [0.0, 0.0]),
        ((-1, -1), [1.0, 0.0]),
        ((0, -1), [0.0, 0.0]),
        ((0, -1), [1.0, 0.0])
    ]

    # Test the result
    @test combined_result == expected_result

    # Test with a different number of nodes (1 node)
    nodes2 = [[0.0, 0.0]]
    nodes2 = reshape(nodes2, :, 1)

    combined_result2 = wavenumber_func.combine_wavenumber_with_all_nodes(wavenumbers, nodes2)
    println(combined_result2)
    expected_result2 = [
        ((-1, -1), [0.0, 0.0]),
        ((0, -1), [0.0, 0.0]),
    ]
    @test combined_result2 == expected_result2

end
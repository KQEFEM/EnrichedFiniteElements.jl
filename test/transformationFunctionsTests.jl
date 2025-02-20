module transformationFunctionsTests

using Test
using EnrichedFiniteElements
const transformations = EnrichedFiniteElements.transformationFunctions

@testset "correct_triangle_orientation Tests" begin
    @testset "Incorrect Orientation - Correction Applied" begin
        # Test case: Incorrectly oriented triangle (clockwise in 2D for negative cross product)
        triangle_nodes_incorrect = [
            0.0 0.0;
            0.5 1.0; # Point 2 is "higher" and to the right, creating clockwise order from 1->2->3
            1.0 0.0;
        ]
        connect_triangle_incorrect = [1, 2, 3]
        expected_triangle_nodes_corrected = [
            0.0 0.0;
            1.0 0.0; # Point 2 and 3 swapped to get counter-clockwise order
            0.5 1.0;
        ]
        expected_connect_triangle_corrected = [1, 3, 2] # Corresponding connectivity update

        updated_triangle_nodes, updated_connect_triangle = transformations.correct_triangle_orientation!(deepcopy(triangle_nodes_incorrect), deepcopy(connect_triangle_incorrect)) # deepcopy to avoid modifying original test data

        @test updated_triangle_nodes ≈ expected_triangle_nodes_corrected # Use ≈ for floating-point comparisons
        @test updated_connect_triangle == expected_connect_triangle_corrected
    end

    @testset "Correct Orientation - No Correction" begin
        # Test case: Correctly oriented triangle (counter-clockwise in 2D for non-negative cross product)
        triangle_nodes_correct = [
            0.0 0.0;
            1.0 0.0; # Point 2 is to the right, point 3 is "above", counter-clockwise order
            0.5 1.0;
        ]
        connect_triangle_correct = [1, 2, 3]

        updated_triangle_nodes, updated_connect_triangle = transformations.correct_triangle_orientation!(deepcopy(triangle_nodes_correct), deepcopy(connect_triangle_correct))

        @test updated_triangle_nodes == triangle_nodes_correct # Should be identical if no correction
        @test updated_connect_triangle == connect_triangle_correct
    end

end

end # module transformationFunctionsTests
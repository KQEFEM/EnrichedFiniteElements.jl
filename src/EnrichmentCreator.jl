module EnrichmentCreator
using LinearAlgebra # Needed for norm
using Base.Iterators

"""
    create_wavenumbers(x_enrichments::Real, y_enrichments::Real)

    Creates a matrix of wavenumber combinations.

    Args:
        x_enrichments: The number of enrichments in the x-direction.
        y_enrichments: The number of enrichments in the y-direction.

    Returns:
        A matrix containing all possible combinations of wavenumbers within the specified range.
"""
function create_wavenumbers(x_enrichments::Real, y_enrichments::Real)
    wavenumbers = collect(Iterators.product(-x_enrichments:x_enrichments, -y_enrichments:y_enrichments))
    return reshape(wavenumbers,:,1)
end

"""
    wavenumber_creation(wavenumbers_ansatz::Matrix{<:Real}, wavenumbers_test::Matrix{<:Real}, Number_of_Frequencies_per_Wavenumber::Int64 = 2)

    Creates a matrix of wavenumber-frequency combinations for ansatz and test functions.

    Args:
        wavenumbers_ansatz: Matrix of wavenumbers for the ansatz functions.
        wavenumbers_test: Matrix of wavenumbers for the test functions.
        Number_of_Frequencies_per_Wavenumber: The number of frequencies to consider per wavenumber. 
                                        Defaults to 2 (positive and negative frequencies).

    Returns:
        A tuple containing:
            - wavenumber_frequency_matrix_ansatz: Matrix of wavenumber-frequency combinations for ansatz functions.
            - wavenumber_frequency_matrix_test: Matrix of wavenumber-frequency combinations for test functions.
"""
function wavenumber_creation(wavenumbers_ansatz::Matrix{<:Real}, wavenumbers_test::Matrix{<:Real}, Number_of_Frequencies_per_Wavenumber::Int64 = 2)
    num_rows = size(wavenumbers_ansatz, 1) * Number_of_Frequencies_per_Wavenumber
    wavenumber_frequency_matrix_ansatz = zeros(Float64, num_rows, 3) 
    wavenumber_frequency_matrix_test = zeros(Float64, num_rows, 3) 

    for ii = 1:size(wavenumbers_ansatz, 1)
        freq_ansatz = norm(wavenumbers_ansatz[ii, :])  
        freq_test = norm(wavenumbers_test[ii, :])      

        if Number_of_Frequencies_per_Wavenumber == 1
            wavenumber_frequency_matrix_ansatz[ii, :] = [wavenumbers_ansatz[ii, 1], wavenumbers_ansatz[ii, 2], 0.0] 
            wavenumber_frequency_matrix_test[ii, :] = [wavenumbers_test[ii, 1], wavenumbers_test[ii, 2], 0.0] 
        else
            wavenumber_frequency_matrix_ansatz[(ii - 1) * 2 + 1, :] = [wavenumbers_ansatz[ii, 1], wavenumbers_ansatz[ii, 2], -freq_ansatz]
            wavenumber_frequency_matrix_ansatz[(ii - 1) * 2 + 2, :] = [wavenumbers_ansatz[ii, 1], wavenumbers_ansatz[ii, 2], freq_ansatz]

            wavenumber_frequency_matrix_test[(ii - 1) * 2 + 1, :] = [wavenumbers_test[ii, 1], wavenumbers_test[ii, 2], -freq_test]
            wavenumber_frequency_matrix_test[(ii - 1) * 2 + 2, :] = [wavenumbers_test[ii, 1], wavenumbers_test[ii, 2], freq_test]
        end
    end

    wavenumber_frequency_matrix_ansatz = unique(wavenumber_frequency_matrix_ansatz, dims=1)
    wavenumber_frequency_matrix_test = unique(wavenumber_frequency_matrix_test, dims=1)

    return wavenumber_frequency_matrix_ansatz, wavenumber_frequency_matrix_test
end

function combine_wavenumber_with_all_nodes(matrix_1::Matrix{T}, matrix_2::Matrix{U}) where {T,U}
    """
    Creates a matrix of all possible combinations of rows from two input matrices.

    Args:
        matrix_1: The first input matrix.
        matrix_2: The second input matrix.

    Returns:
        A matrix containing all possible combinations of rows from the input matrices.

    Raises:
        ArgumentError: If either `matrix_1` or `matrix_2` is empty.
    """
    if isempty(matrix_1) || isempty(matrix_2)
        throw(ArgumentError("Both input matrices must have at least one row."))
    end

    # Create an iterator of all combinations
    combinations = product(eachrow(matrix_1), eachrow(matrix_2))

    # Convert the iterator to an array of tuples
    result = collect(combinations)

    return result
end

end 
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
#@docs EnrichedFiniteElements.EnrichmentCreator.create_wavenumbers
function create_wavenumbers(x_enrichments::Real, y_enrichments::Real)
    wavenumbers = reshape(
        collect(
            Iterators.product(-x_enrichments:x_enrichments, -y_enrichments:y_enrichments),
        ),
        :,
        1,
    )
    return hcat(first.(wavenumbers), last.(wavenumbers))

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
#@docs EnrichedFiniteElements.EnrichmentCreator.wavenumber_creation

function wavenumber_creation(
    wavenumbers_ansatz::Matrix{<:Real},
    wavenumbers_test::Matrix{<:Real},
    Number_of_Frequencies_per_Wavenumber::Int64 = 2,
)
    num_rows = size(wavenumbers_ansatz, 1) * Number_of_Frequencies_per_Wavenumber
    wavenumber_frequency_matrix_ansatz = zeros(Float64, num_rows, 3)
    wavenumber_frequency_matrix_test = zeros(Float64, num_rows, 3)

    for ii = 1:size(wavenumbers_ansatz, 1)
        freq_ansatz = norm(wavenumbers_ansatz[ii, :])
        freq_test = norm(wavenumbers_test[ii, :])

        if Number_of_Frequencies_per_Wavenumber == 1
            wavenumber_frequency_matrix_ansatz[ii, :] =
                [wavenumbers_ansatz[ii, 1], wavenumbers_ansatz[ii, 2], 0.0]
            wavenumber_frequency_matrix_test[ii, :] =
                [wavenumbers_test[ii, 1], wavenumbers_test[ii, 2], 0.0]
        else
            wavenumber_frequency_matrix_ansatz[(ii-1)*2+1, :] =
                [wavenumbers_ansatz[ii, 1], wavenumbers_ansatz[ii, 2], -freq_ansatz]
            wavenumber_frequency_matrix_ansatz[(ii-1)*2+2, :] =
                [wavenumbers_ansatz[ii, 1], wavenumbers_ansatz[ii, 2], freq_ansatz]

            wavenumber_frequency_matrix_test[(ii-1)*2+1, :] =
                [wavenumbers_test[ii, 1], wavenumbers_test[ii, 2], -freq_test]
            wavenumber_frequency_matrix_test[(ii-1)*2+2, :] =
                [wavenumbers_test[ii, 1], wavenumbers_test[ii, 2], freq_test]
        end
    end
    # **Normalize Zero Values to Ensure -0.0 and 0.0 Are Treated as the Same**
    function normalize_zeros(matrix)
        return map(x -> iszero(x) ? 0.0 : x, matrix)
    end
    wavenumber_frequency_matrix_ansatz = normalize_zeros(wavenumber_frequency_matrix_ansatz)
    wavenumber_frequency_matrix_test = normalize_zeros(wavenumber_frequency_matrix_test)

    wavenumber_frequency_matrix_ansatz =
    sortslices(unique(wavenumber_frequency_matrix_ansatz, dims = 1), dims=1, by=x -> x[1])
    wavenumber_frequency_matrix_test = sortslices(unique(wavenumber_frequency_matrix_test, dims = 1), dims=1, by=x -> x[1])

    return wavenumber_frequency_matrix_ansatz, wavenumber_frequency_matrix_test
end

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
#@docs EnrichedFiniteElements.EnrichmentCreator.combine_wavenumber_with_all_nodes

function combine_wavenumber_with_all_nodes(
    matrix_1::Matrix{T},
    matrix_2::Matrix{U},
) where {T,U}

    if isempty(matrix_1) || isempty(matrix_2)
        throw(ArgumentError("Both input matrices must have at least one row."))
    end

    # Create an iterator of all combinations
    combinations = product(eachrow(matrix_1), eachrow(matrix_2))

    # Convert the iterator to an array of tuples
    result = collect(combinations)
    result = reshape(result, :, 1)

    return result
end
""" 
Creates the pairs of indices for a connectivity matrix for test and ansatz
Notice that the ordering is slightly differnet in the julia version compared to matlab. The wavenumber idx goes [1,1],[2,1] compared to matlab [1,1],[1,2]
"""
function generate_pairs(vector_1::AbstractVector{T}, vector_2::AbstractVector{T}) where {T}
    if isempty(vector_1) || isempty(vector_2)
        throw(ArgumentError("The input vectors must not be empty."))
    end

    # Use broadcasting to generate all pairs
    pairs = [(x, y) for x in vector_1, y in vector_2]
    pairs = reshape(pairs, :, 1)
    return pairs
end

end

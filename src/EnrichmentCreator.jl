module EnrichmentCreator
using LinearAlgebra # Needed for norm
using Base.Iterators

function create_wavenumbers(x_enrichments::Real, y_enrichments::Real)
    """ Creates the wavenumber matrix to use """
    wavenumbers = collect(Iterators.product(-x_enrichments:x_enrichments, -y_enrichments:y_enrichments))
    return reshape(wavenumbers,:,1)
    
end 

function wavenumber_creation(wavenumbers_ansatz::Matrix{<:Real}, wavenumbers_test::Matrix{<:Real}, Number_of_Frequencies_per_Wavenumber::Int64 = 2)
    num_rows = size(wavenumbers_ansatz, 1) * Number_of_Frequencies_per_Wavenumber
    wavenumber_frequency_matrix_ansatz = zeros(Float64, num_rows, 3) # Specify Float64
    wavenumber_frequency_matrix_test = zeros(Float64, num_rows, 3) # Specify Float64

    for ii = 1:size(wavenumbers_ansatz, 1)
        freq_ansatz = norm(wavenumbers_ansatz[ii, :])  # Calculate norm once
        freq_test = norm(wavenumbers_test[ii, :])      # Calculate norm once

        if Number_of_Frequencies_per_Wavenumber == 1
            wavenumber_frequency_matrix_ansatz[ii, :] = [wavenumbers_ansatz[ii, 1], wavenumbers_ansatz[ii, 2], 0.0] # Correct frequency
            wavenumber_frequency_matrix_test[ii, :] = [wavenumbers_test[ii, 1], wavenumbers_test[ii, 2], 0.0] # Correct frequency
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
    # Create an iterator of all combinations
    combinations = product(eachrow(matrix_1), eachrow(matrix_2))

    # Convert the iterator to an array of tuples
    result = collect(combinations)

    return result
end

end 
function combine_wavenumber_with_all_nodes(wavenumbers, nodes)
    n_wavenumbers = size(wavenumbers, 1)
    n_nodes = size(nodes, 1)

    result = Array{Tuple{Vector{eltype(wavenumbers)}, Vector{eltype(nodes)}}}(undef, n_wavenumbers * n_nodes)

    for k_iter in 1:n_wavenumbers
        wavenumber = wavenumbers[k_iter, :]  # Current wavenumber row
        for node_iter in 1:n_nodes
          result[(k_iter-1)*n_nodes + node_iter] = (wavenumber, nodes[node_iter, :])
        end
    end

    return result
end


# Example usage:
n_wavenumbers = 3
n_nodes = 5
wavenumbers = rand(n_wavenumbers, 3)
nodes = rand(n_nodes, 2)
wavenumbers = [-1,0,1]; nodes =[-1,0,1];
combined = combine_wavenumber_with_all_nodes(wavenumbers, nodes)
c = combine_wavenumber_with_all_nodes(combined, combined)

# Print or process the combinations as needed
for (wavenumber, node) in combined
    println("Wavenumber: $wavenumber, Node: $node")
end

println("Number of combinations: ", length(combined)) #Should be n_wavenumbers * n_nodes

# Accessing a specific combination:
wavenumber_2_node_3 = combined[(2-1)*n_nodes + 3] #Wavenumber 2 with node 3
println("Wavenumber 2, Node 3: Wavenumber $(wavenumber_2_node_3[1]), Node $(wavenumber_2_node_3[2])")

wavenumber_idx = 2
node_idx = 4
wavenumber_node_combination = combined[(wavenumber_idx-1)*n_nodes + node_idx]
println("Wavenumber $wavenumber_idx, Node $node_idx: Wavenumber $(wavenumber_node_combination[1]), Node $(wavenumber_node_combination[2])")


wavenumbers = collect(Iterators.product(-2:1, -1:1))
wavenumbers = reshape(wavenumbers,:,1)
combined = combine_wavenumber_with_all_nodes(wavenumbers, nodes)


# Initialize variables
using LinearAlgebra  # Import the LinearAlgebra module

# Initialize variables
NO_Frequencies = 0  # Set this to the actual value
Number_of_Frequencies_per_Wavenumber = 2  # Set this to the actual value

# Define your wavenumbers (each row is a 2D vector)
wavenumber_ansatz = [
    -1 -1;
    0 -1;
    1 -1;
    -1  0;
    0  0;
    1  0;
    -1  1;
    0  1;
    1  1;
]

wavenumber_test = wavenumber_ansatz  # Assuming wavenumber_test is the same for simplicity

# Initialize empty lists for Freq_Ansatz and Freq_Test
Freq_Ansatz = [zeros(2) for _ in 1:size(wavenumber_ansatz, 1)]
Freq_Test = [zeros(2) for _ in 1:size(wavenumber_ansatz, 1)]

# Initialize a list to store pairs of (wavenumber, frequency)
wavenumber_frequency_pairs_ansatz = []
wavenumber_frequency_pairs_test = []

# Loop through the wavenumbers
for ii = 1:size(wavenumber_ansatz, 1)
    if NO_Frequencies != 1
        if Number_of_Frequencies_per_Wavenumber == 2
            # Calculate the frequencies
            freq_ansatz_1 = -1 * norm(wavenumber_ansatz[ii, :])
            freq_ansatz_2 = 1 * norm(wavenumber_ansatz[ii, :])

            freq_test_1 = -1 * norm(wavenumber_test[ii, :])
            freq_test_2 = 1 * norm(wavenumber_test[ii, :])
            
            # Add the flattened (wavenumber, frequency) pairs to the lists
            push!(wavenumber_frequency_pairs_ansatz, (wavenumber_ansatz[ii, 1], wavenumber_ansatz[ii, 2], freq_ansatz_1))
            push!(wavenumber_frequency_pairs_ansatz, (wavenumber_ansatz[ii, 1], wavenumber_ansatz[ii, 2], freq_ansatz_2))
            
            push!(wavenumber_frequency_pairs_test, (wavenumber_test[ii, 1], wavenumber_test[ii, 2], freq_test_1))
            push!(wavenumber_frequency_pairs_test, (wavenumber_test[ii, 1], wavenumber_test[ii, 2], freq_test_2))

        end
    else
        push!(wavenumber_frequency_pairs_test, (wavenumber_test[ii, 1], wavenumber_test[ii, 2], 0))
        push!(wavenumber_frequency_pairs_ansatz, (wavenumber_ansatz[ii, 1], wavenumber_ansatz[ii, 2], 0))
    end
end

# Optional: print or inspect the result
println("Flattened Wavenumber-Frequency pairs for Ansatz:")
println(wavenumber_frequency_pairs_ansatz)

println("Flattened Wavenumber-Frequency pairs for Test:")
println(wavenumber_frequency_pairs_test)
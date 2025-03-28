using IterTools

function combine_wavenumber_with_all_nodes(
    matrix_1::Matrix{T},
    matrix_2::Matrix{U},
) where {T,U}
    if isempty(matrix_1) || isempty(matrix_2)
        throw(ArgumentError("Both input matrices must have at least one row."))
    end

    # Create an iterator of all combinations of rows
    combinations = product(eachrow(matrix_1), eachrow(matrix_2))

    # Convert the iterator to an array of tuples
    result = collect(combinations)

    return result
end


# Example usage:
n_wavenumbers = 3
n_nodes = 5
wavenumbers = rand(n_wavenumbers, 3)
nodes = rand(n_nodes, 2)
wavenumbers = [-1, 0, 1];
nodes = [-1, 0, 1];
combined = combine_wavenumber_with_all_nodes(wavenumbers, nodes)
c = combine_wavenumber_with_all_nodes(combined, combined)

# Print or process the combinations as needed
for (wavenumber, node) in combined
    println("Wavenumber: $wavenumber, Node: $node")
end

println("Number of combinations: ", length(combined)) #Should be n_wavenumbers * n_nodes

# Accessing a specific combination:
wavenumber_2_node_3 = combined[(2-1)*n_nodes+3] #Wavenumber 2 with node 3
println(
    "Wavenumber 2, Node 3: Wavenumber $(wavenumber_2_node_3[1]), Node $(wavenumber_2_node_3[2])",
)

wavenumber_idx = 2
node_idx = 4
wavenumber_node_combination = combined[(wavenumber_idx-1)*n_nodes+node_idx]
println(
    "Wavenumber $wavenumber_idx, Node $node_idx: Wavenumber $(wavenumber_node_combination[1]), Node $(wavenumber_node_combination[2])",
)


wavenumbers = collect(Iterators.product(-2:1, -1:1))
wavenumbers = reshape(wavenumbers, :, 1)
combined = combine_wavenumber_with_all_nodes(wavenumbers, nodes)



# Initialize variables
using LinearAlgebra  # Import the LinearAlgebra module

# Initialize variables
NO_Frequencies = 0  # Set this to the actual value
Number_of_Frequencies_per_Wavenumber = 2  # Set this to the actual value

# Define your wavenumbers (each row is a 2D vector)
wavenumbers_ansatz = [
    -1 -1
    0 -1
    1 -1
    -1 0
    0 0
    1 0
    -1 1
    0 1
    1 1
]
wavenumbers_ansatz = Union{Int64,Float64}[1 1; 1.0 -1]
wavenumbers_test = wavenumbers_ansatz
num_rows = size(wavenumbers_ansatz, 1) * Number_of_Frequencies_per_Wavenumber
wavenumber_frequency_matrix_ansatz = zeros(Real, num_rows, 3)
wavenumber_frequency_matrix_test = zeros(Real, num_rows, 3)

# Loop through the wavenumbers and calculate frequencies directly
num_rows = size(wavenumbers_ansatz, 1) * Number_of_Frequencies_per_Wavenumber
wavenumber_frequency_matrix_ansatz = zeros(Real, num_rows, 3)
wavenumber_frequency_matrix_test = zeros(Real, num_rows, 3)

# Loop through the wavenumbers and calculate frequencies directly
for ii = 1:size(wavenumbers_ansatz, 1)

    # Calculate the frequencies for Ansatz and Test
    freq_ansatz_1 = -1 * norm(wavenumbers_ansatz[ii, :])
    freq_ansatz_2 = 1 * norm(wavenumbers_ansatz[ii, :])

    freq_test_1 = -1 * norm(wavenumbers_test[ii, :])
    freq_test_2 = 1 * norm(wavenumbers_test[ii, :])

    # Add the (wavenumber, frequency) pairs to the matrix
    wavenumber_frequency_matrix_ansatz[(ii-1)*2+1, :] .=
        [wavenumbers_ansatz[ii, 1], wavenumbers_ansatz[ii, 2], freq_ansatz_1]
    wavenumber_frequency_matrix_ansatz[(ii-1)*2+2, :] .=
        [wavenumbers_ansatz[ii, 1], wavenumbers_ansatz[ii, 2], freq_ansatz_2]

    wavenumber_frequency_matrix_test[(ii-1)*2+1, :] .=
        [wavenumbers_test[ii, 1], wavenumbers_test[ii, 2], freq_test_1]
    wavenumber_frequency_matrix_test[(ii-1)*2+2, :] .=
        [wavenumbers_test[ii, 1], wavenumbers_test[ii, 2], freq_test_2]



end

wavenumber_frequency_matrix_ansatz = unique(wavenumber_frequency_matrix_ansatz, dims = 1)
wavenumber_frequency_matrix_test = unique(wavenumber_frequency_matrix_test, dims = 1)
wavenumber_frequency_matrix_ansatz = unique(wavenumber_frequency_matrix_ansatz, dims = 1)
wavenumber_frequency_matrix_test = unique(wavenumber_frequency_matrix_test, dims = 1)

# Optional: print or inspect the result
println("Flattened Wavenumber-Frequency pairs for Ansatz:")
println(wavenumber_frequency_matrix_ansatz)

println("Flattened Wavenumber-Frequency pairs for Test:")
println(wavenumber_frequency_matrix_test)

a = [1.1 2.3 3; 1 1.1 2; 1.1 2.3 3; 2.3 3 1]
unique(a, dims = 1)

using Revise
using EnrichedFiniteElements


const wave_func = EnrichedFiniteElements.EnrichmentCreator
wave_func.create_wavenumbers(1, 1)


matrix_1 = rand(10, 3)  # Example 10x3 matrix
matrix_2 = rand(5, 2)   # Example 5x2 matrix

result = wave_func.combine_wavenumber_with_all_nodes(matrix_1, matrix_2)

# This will be for the looping

domain = ((0, 1), (0, 1))
num_nodes = 4
const mesh_create = EnrichedFiniteElements.MeshCreation

# Create the mesh
mesh = mesh_create.rectangle_domain(domain)

nodes = mesh.nodes
connectivity = mesh.connectivity
boundary_index = mesh.boundary_idx
boundary_edges = mesh.boundary_edges
nodes_idx = collect(1:size(nodes, 1))
wavenumber_idx_ansatz = [1 2 3 4 5 6 7]#reshape(collect(1:size(matrix_1, 1)),(1,10))
wavenumber_idx_test = wavenumber_idx_ansatz
result_wavenumber_idx =
    wave_func.combine_wavenumber_with_all_nodes(wavenumber_idx_ansatz, wavenumber_idx_test)

#! This can be used to vectorise the computation as we just loop through this list
wavenumber_connect_idx =
    a = wave_func.combine_wavenumber_with_all_nodes(result_wavenumber_idx, connectivity)


wavenumber_idx_ansatz = reshape(wavenumber_idx_ansatz, 1, :)

#! This is the code to use for the matrix creation
using Revise
using EnrichedFiniteElements
domain = ((0, 1), (0, 1))
num_nodes = 4
const mesh_create = EnrichedFiniteElements.MeshCreation
const wave_func = EnrichedFiniteElements.EnrichmentCreator

# Create the mesh
mesh = mesh_create.rectangle_domain(domain)

nodes = mesh.nodes
connectivity = mesh.connectivity
boundary_index = mesh.boundary_idx
boundary_edges = mesh.boundary_edges

vector = [1, 2, 3, 4, 5, 6, 7]
vector_2 = [10, 20, 30, 40, 50, 60, 70]
vector = [1, 2, 3]
vector_2 = [60, 70]
all_pairs = wave_func.generate_pairs(vector, vector_2)
println(all_pairs)
result = wave_func.combine_wavenumber_with_all_nodes(all_pairs, connectivity)






domain = ((0, 1), (0, 1))
num_nodes = 2
const mesh_create = EnrichedFiniteElements.MeshCreation

# Create the mesh
mesh = mesh_create.rectangle_domain(domain)
connectivity = [1 2 3; 4 5 6]  # 2x3 matrix with sequential values
vector = [1, 2, 3]
vector_2 = [60, 70]

all_pairs = wave_func.generate_pairs(vector, vector_2)
println(all_pairs)
result = wave_func.combine_wavenumber_with_all_nodes(all_pairs, connectivity)

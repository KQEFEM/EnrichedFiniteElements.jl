using Revise
using LinearAlgebra  # Import the LinearAlgebra module
using SparseArrays
using EnrichedFiniteElements

const Transformations = EnrichedFiniteElements.TransformationFunctions
const integrator = EnrichedFiniteElements.Operators

const matrix_comp = EnrichedFiniteElements.MatrixCreation

const mesh_create = EnrichedFiniteElements.MeshCreation
const wave_func = EnrichedFiniteElements.EnrichmentCreator

#! Set up for the matrix creation
domain = ((0, 1), (0, 1))
num_nodes = 10

# Create the mesh
mesh = mesh_create.rectangle_domain(domain)

nodes = mesh.nodes
connectivity = mesh.connectivity
boundary_index = mesh.boundary_idx
boundary_edges = mesh.boundary_edges

wave_x = 1
wave_y = 1
ansatz_wave = wave_func.create_wavenumbers(wave_x, wave_y)
test_ansatz = wave_func.create_wavenumbers(wave_x, wave_y)
wavenumbers_ansatz, wavenumbers_test =
    wave_func.wavenumber_creation(ansatz_wave, test_ansatz)

idx_wave_ansatz = collect(1:size(wavenumbers_ansatz, 1))
idx_wave_test = collect(1:size(wavenumbers_test, 1))

all_pairs = wave_func.generate_pairs(idx_wave_ansatz, idx_wave_test)
println(all_pairs)
idx_connectivity = collect(1:size(connectivity, 1))
result = wave_func.combine_wavenumber_with_all_nodes(all_pairs, connectivity)

wave_node_pairs = result
dt = 0.1
using Revise
mass_cell, _ = matrix_comp.compute_sparse_matrix(
    all_pairs,
    nodes,
    wave_node_pairs,
    wavenumbers_ansatz,
    wavenumbers_test,
    integrator,
    dt,
    t_jump = dt,
    t0 = 0.0,
    mass_bool = true,
)
mass_cell

# Assemble the final large matrix without loops
num_rows = size(cell_sparse_zero_array_2, 1)
num_cols = size(cell_sparse_zero_array_2, 2)
final_matrix = reduce(
    vcat,
    [reduce(hcat, [cell_sparse_zero_array_2[i, j] for j = 1:num_cols]) for i = 1:num_rows],
);

exact_mat = [
    0.031250000000000 0.000000000000000 0.000000000000000 0.000000000000000 0.007812500000000 0.000000000000000 0.000000000000000 0.007812500000000 0.000000000000000 0.015625000000000 0.000000000000000 0.000000000000000
    0.000000000000000 0.023437500000000 0.000000000000000 0.000000000000000 0.005859375000000 0.005859375000000 0.000000000000000 0.000000000000000 0.000000000000000 0.000000000000000 0.000000000000000 0.011718750000000
    0.000000000000000 0.000000000000000 0.029513875000000 0.000000000000000 0.000000000000000 0.007335062500000 0.007421875000000 0.000000000000000 0.000000000000000 0.000000000000000 0.014756937500000 0.000000000000000
    0.000000000000000 0.000000000000000 0.000000000000000 0.024479166666667 0.000000000000000 0.000000000000000 0.006119791666667 0.006119791666667 0.012239583333333 0.000000000000000 0.000000000000000 0.000000000000000
    0.007812500000000 0.005859375000000 0.000000000000000 0.000000000000000 0.037109375000000 0.000000000000000 0.000000000000000 0.000000000000000 0.000000000000000 0.012695312500000 0.000000000000000 0.010742187500000
    0.000000000000000 0.005859375000000 0.007335062500000 0.000000000000000 0.000000000000000 0.036176195312500 0.000000000000000 0.000000000000000 0.000000000000000 0.000000000000000 0.012228722656250 0.010753035156250
    0.000000000000000 0.000000000000000 0.007421875000000 0.006119791666667 0.000000000000000 0.000000000000000 0.036827265104167 0.000000000000000 0.010991757552083 0.000000000000000 0.012293840885417 0.000000000000000
    0.007812500000000 0.000000000000000 0.000000000000000 0.006119791666667 0.000000000000000 0.000000000000000 0.000000000000000 0.037369791666667 0.010872395833333 0.012565104166667 0.000000000000000 0.000000000000000
    0.000000000000000 0.000000000000000 0.000000000000000 0.012239583333333 0.000000000000000 0.000000000000000 0.010991757552083 0.010872395833333 0.053081614583333 0.009429258072917 0.009548619791667 0.000000000000000
    0.015625000000000 0.000000000000000 0.000000000000000 0.000000000000000 0.012695312500000 0.000000000000000 0.000000000000000 0.012565104166667 0.009429258072917 0.069704872916667 0.009592019791667 0.009798178385417
    0.000000000000000 0.000000000000000 0.014756937500000 0.000000000000000 0.000000000000000 0.012228722656250 0.012293840885417 0.000000000000000 0.009548619791667 0.009592019791667 0.068229166666667 0.009809026041667
    0.000000000000000 0.011718750000000 0.000000000000000 0.000000000000000 0.010742187500000 0.010753035156250 0.000000000000000 0.000000000000000 0.000000000000000 0.009798178385417 0.009809026041667 0.052821177083333
];
exact_mat = [
    0.030877602472284-0.003877008570769im 0.000000000000000 0.000000000000000 0.000000000000000 0.007775970724166-0.000583979482525im 0.000000000000000 0.000000000000000 0.007486048785433-0.002112407182345im 0.000000000000000 0.015245131004788-0.003089998063468im 0.000000000000000 0.000000000000000
    0.000000000000000 0.023224095283954-0.002474700051558im 0.000000000000000 0.000000000000000 0.005843947279790-0.000328969785424im 0.005644388243485-0.001478931801369im 0.000000000000000 0.000000000000000 0.000000000000000 0.000000000000000 0.000000000000000 0.011523333568294-0.001888668887620im
    0.000000000000000 0.000000000000000 0.018749118824709-0.022640302346737im 0.000000000000000 0.000000000000000 0.005448546640172-0.004863509288037im 0.004436764991800-0.005934077088816im 0.000000000000000 0.000000000000000 0.000000000000000 0.010154738647367-0.010620645932694im 0.000000000000000
    0.000000000000000 0.000000000000000 0.000000000000000 0.015327008867588-0.018974744895991im 0.000000000000000 0.000000000000000 0.003598961356286-0.004940977585399im 0.004495262533220-0.004114807165160im 0.008207089397093-0.009019825982433im 0.000000000000000 0.000000000000000 0.000000000000000
    0.007775970724166-0.000583979482525im 0.005843947279790-0.000328969785424im 0.000000000000000 0.000000000000000 0.036899953229800-0.003100947021406im 0.000000000000000 0.000000000000000 0.000000000000000 0.000000000000000 0.012470658492625-0.002161604280700im 0.000000000000000 0.010600531331030-0.001565045746753im
    0.000000000000000 0.005644388243485-0.001478931801369im 0.005448546640172-0.004863509288037im 0.000000000000000 0.000000000000000 0.031424206939443-0.017232543150414im 0.000000000000000 0.000000000000000 0.000000000000000 0.000000000000000 0.010045566118043-0.006871466278405im 0.009977411412215-0.003882939027479im
    0.000000000000000 0.000000000000000 0.004436764991800-0.005934077088816im 0.003598961356286-0.004940977585399im 0.000000000000000 0.000000000000000 0.022341021818855-0.029181810189094im 0.000000000000000 0.007228733521186-0.008244105175587im 0.000000000000000 0.008234898825151-0.009080987790504im 0.000000000000000
    0.007486048785433-0.002112407182345im 0.000000000000000 0.000000000000000 0.004495262533220-0.004114807165160im 0.000000000000000 0.000000000000000 0.000000000000000 0.032529214730831-0.017691336840626im 0.008763060958276-0.006359906505235im 0.011503568856235-0.004910938483760im 0.000000000000000 0.000000000000000
    0.000000000000000 0.000000000000000 0.000000000000000 0.008207089397093-0.009019825982433im 0.000000000000000 0.000000000000000 0.007228733521186-0.008244105175587im 0.008763060958276-0.006359906505235im 0.039958205695791-0.034536070242700im 0.008038388824307-0.004892627400856im 0.007407464468733-0.005975080448259im 0.000000000000000
    0.015245131004788-0.003089998063468im 0.000000000000000 0.000000000000000 0.000000000000000 0.012470658492625-0.002161604280700im 0.000000000000000 0.000000000000000 0.011503568856235-0.004910938483760im 0.008038388824307-0.004892627400856im 0.065024308941807-0.023775853493687im 0.008371025634916-0.004629748403171im 0.009245335433152-0.003136767827467im
    0.000000000000000 0.000000000000000 0.010154738647367-0.010620645932694im 0.000000000000000 0.000000000000000 0.010045566118043-0.006871466278405im 0.008234898825151-0.009080987790504im 0.000000000000000 0.007407464468733-0.005975080448259im 0.008371025634916-0.004629748403171im 0.053641165033880-0.041414029649822im 0.008779916473549-0.004322366957603im
    0.000000000000000 0.011523333568294-0.001888668887620im 0.000000000000000 0.000000000000000 0.010600531331030-0.001565045746753im 0.009977411412215-0.003882939027479im 0.000000000000000 0.000000000000000 0.000000000000000 0.009245335433152-0.003136767827467im 0.008779916473549-0.004322366957603im 0.050421595925460-0.014777582537093im
]

Mass_matlab = [
    0.0312 0 0 0 0.0078 0 0 0.0078 0 0.0156 0 0
    0 0.0234 0 0 0.0059 0.0059 0 0 0 0 0 0.0117
    0 0 0.0295 0 0 0.0073 0.0074 0 0 0 0.0148 0
    0 0 0 0.0245 0 0 0.0061 0.0061 0.0122 0 0 0
    0.0078 0.0059 0 0 0.0371 0 0 0 0 0.0127 0 0.0107
    0 0.0059 0.0073 0 0 0.0362 0 0 0 0 0.0122 0.0108
    0 0 0.0074 0.0061 0 0 0.0368 0 0.0110 0 0.0123 0
    0.0078 0 0 0.0061 0 0 0 0.0374 0.0109 0.0126 0 0
    0 0 0 0.0122 0 0 0.0110 0.0109 0.0531 0.0094 0.0095 0
    0.0156 0 0 0 0.0127 0 0 0.0126 0.0094 0.0697 0.0096 0.0098
    0 0 0.0148 0 0 0.0122 0.0123 0 0.0095 0.0096 0.0682 0.0098
    0 0.0117 0 0 0.0107 0.0108 0 0 0 0.0098 0.0098 0.0528
]

println(norm(cell_sparse_zero_array_2[end-1, end] - conj(exact_mat)))

using DelimitedFiles
using SparseArrays

# Load the data from the text file
data = readdlm("test/testdata/MassMatrixEnriched_enriched.txt")

# Extract row indices, column indices, and values
rows = Int.(data[:, 1])
cols = Int.(data[:, 2])
vals = complex.(data[:, 3], data[:, 4])  # Combine real and imaginary parts

# Reconstruct the sparse matrix
Mass_matlab = sparse(rows, cols, vals);

println(norm((final_matrix - conj(Mass_matlab))))
println(norm((final_matrix - Mass_matlab)))





# Find the index of the row containing exactly [0.0, 0.0]
row_index = findfirst(row -> all(row .== [0.0, 0.0]), eachrow(wavenumbers_ansatz[:, 1:2]))

# Print the result (1-based index)
println("Row location: ", row_index)

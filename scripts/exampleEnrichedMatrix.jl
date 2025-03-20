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
    wave_func.wavenumber_creation(ansatz_wave, test_ansatz, 2)

idx_wave_ansatz = collect(1:size(wavenumbers_ansatz, 1))
idx_wave_test = collect(1:size(wavenumbers_test, 1))

all_pairs = wave_func.generate_pairs(idx_wave_ansatz, idx_wave_test)
println(all_pairs)
idx_connectivity = collect(1:size(connectivity, 1))
result = wave_func.combine_wavenumber_with_all_nodes(all_pairs, connectivity)

wave_node_pairs = result
dt = 0.1
cell_sparse_zero_array_2 = matrix_comp.compute_sparse_mass_matrix(
    all_pairs,
    nodes,
    wave_node_pairs,
    wavenumbers_ansatz,
    wavenumbers_test,
    integrator,
    dt,
    dt,
    0.0,
)

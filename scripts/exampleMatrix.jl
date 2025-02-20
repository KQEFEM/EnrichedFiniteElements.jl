using Revise
using EnrichedFiniteElements

const mesh_create = EnrichedFiniteElements.MeshCreation
const wave_func = EnrichedFiniteElements.EnrichmentCreator


domain = ((0, 1), (0, 1))
num_nodes = 10

# Create the mesh
mesh = mesh_create.rectangle_domain(domain)

nodes = mesh.nodes
connectivity = mesh.connectivity
boundary_index = mesh.boundary_idx
boundary_edges = mesh.boundary_edges

ansatz_wave = wave_func.create_wavenumbers(1,  1)
test_ansatz = wave_func.create_wavenumbers(1,  1)
wavenumbers = wave_func.wavenumber_creation(
    ansatz_wave,
    test_ansatz,
    2,
)
all_pairs = wave_func.generate_pairs(ansatz_wave, test_ansatz)
println(all_pairs)
result = wave_func.combine_wavenumber_with_all_nodes(all_pairs, connectivity)

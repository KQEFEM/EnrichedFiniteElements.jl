
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

wave_x = 0
wave_y = 0
ansatz_wave = wave_func.create_wavenumbers(wave_x, wave_y)
test_ansatz = wave_func.create_wavenumbers(wave_x, wave_y)
wavenumbers_ansatz, wavenumbers_test = wave_func.wavenumber_creation(
    ansatz_wave,
    test_ansatz,
    Number_of_Frequencies_per_Wavenumber = 2,
    time_enrichment_only = false,
)

idx_wave_ansatz = collect(1:size(wavenumbers_ansatz, 1))
idx_wave_test = collect(1:size(wavenumbers_test, 1))

all_pairs = wave_func.generate_pairs(idx_wave_ansatz, idx_wave_test)
println(all_pairs)
idx_connectivity = collect(1:size(connectivity, 1))
result = wave_func.combine_wavenumber_with_all_nodes(all_pairs, connectivity)

wave_node_pairs = result
dt = 0.1

using Revise
 _, pDtq_cell, vDxeta_cell = matrix_comp.compute_sparse_matrix(
    all_pairs,
    nodes,
    wave_node_pairs,
    wavenumbers_ansatz,
    wavenumbers_test,
    integrator,
    dt,
    convection_bool = true,
)
vDxeta_array = permutedims(vDxeta_cell, (2,1)) #!There is a missing transpose somewhere in the basis operations

pDtq_array = matrix_comp.convert_sparse_cell_to_array(pDtq_cell)
vDxeta_array = matrix_comp.convert_sparse_cell_to_array(vDxeta_cell)
# pDtq_exact = spzeros(size(array, 1), size(array, 2))

using DelimitedFiles
using SparseArrays

# Load the data from the text file
data = readdlm("test/testdata/ConvectionDt_space_time.txt")

# Extract row indices, column indices, and values
rows = Int.(data[:, 1])
cols = Int.(data[:, 2])
vals = complex.(data[:, 3], data[:, 4])  # Combine real and imaginary parts

# Reconstruct the sparse matrix
pDtq_exact = sparse(rows, cols, vals);

println(norm((pDtq_array - transpose(pDtq_exact))))

## Dx Convection 
# Load the data from the text file
data = readdlm("test/testdata/ConvectionDx_FEM.txt")

# Extract row indices, column indices, and values
rows = Int.(data[:, 1])
cols = Int.(data[:, 2])
vals = complex.(data[:, 3], data[:, 4])  # Combine real and imaginary parts

# Reconstruct the sparse matrix
vDxeta_exact = sparse(rows, cols, vals);

println(norm(( vDxeta_array - transpose( vDxeta_exact))))
println(norm(diag( vDxeta_array - ( vDxeta_exact))))

vDxeta_array[1]
vDxeta_exact[1]
vDxeta_array[1,13]
vDxeta_exact[1,13]






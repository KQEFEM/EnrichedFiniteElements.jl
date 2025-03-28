
using Revise
using LinearAlgebra  # Import the LinearAlgebra module
using SparseArrays
using EnrichedFiniteElements

const Transformations = EnrichedFiniteElements.TransformationFunctions
const integrator = EnrichedFiniteElements.Operators

const rhs_comp = EnrichedFiniteElements.LoadTerms

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
rhs_function = rhs_comp.exact_pressure_0_velocity(1.9,2.0,1.0)  # Assign the function to a variable
typeof(rhs_function)
# Call the load_term_integration function
f1, f2x, f2y = rhs_comp.load_term_integration(
    all_pairs,
    nodes,
    wave_node_pairs,
    wavenumbers_ansatz,
    wavenumbers_test;
    dt = dt,
    t0 = 0,
    rhs_function=rhs_function,  # Pass rhs_function as a keyword argument
    kx = 1.0,
    ky = 1.0,
    w = 1.0
)





vDxeta_cell = convection_cell[2]
pDtq_cell = convection_cell[1]
vDyeta_cell = convection_cell[3]
vDxeta_cell = permutedims(vDxeta_cell, (2, 1)) #!There is a missing transpose somewhere in the basis operations
pDtq_cell = permutedims(pDtq_cell, (2, 1)) #!There is a missing transpose somewhere in the basis operations
vDyeta_cell = permutedims(vDyeta_cell, (2, 1))
pDtq_array = matrix_comp.convert_sparse_cell_to_array(pDtq_cell);
vDxeta_array = matrix_comp.convert_sparse_cell_to_array(vDxeta_cell);
vDyeta_array = matrix_comp.convert_sparse_cell_to_array(vDyeta_cell);
C = [pDtq_array vDxeta_array vDyeta_array;
     vDxeta_array pDtq_array pDtq_array;
     vDyeta_array pDtq_array pDtq_array]
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

println(norm((pDtq_array - (pDtq_exact))))
println(norm(imag(pDtq_array - (pDtq_exact))))
println(norm(real(pDtq_array - (pDtq_exact))))
pDtq_array[1]
pDtq_exact[1]
pDtq_array[1, 13]
pDtq_exact[13, 1]
## Dx Convection 
# Load the data from the text file
data = readdlm("test/testdata/ConvectionDx_space_time.txt")

# Extract row indices, column indices, and values
rows = Int.(data[:, 1])
cols = Int.(data[:, 2])
vals = complex.(data[:, 3], data[:, 4])  # Combine real and imaginary parts

# Reconstruct the sparse matrix
vDxeta_exact = (sparse(rows, cols, vals));

println(norm((vDxeta_array - (vDxeta_exact))))
println(norm(real(vDxeta_array - (vDxeta_exact))))
println(norm(imag(vDxeta_array - (vDxeta_exact))))


println(norm(diag(vDxeta_array - (vDxeta_exact))))

vDxeta_array[1]
vDxeta_exact[1]
vDxeta_array[1, 17]
vDxeta_exact[17, 1]





vDxeta_array[1, 17] - vDxeta_exact[17, 1]

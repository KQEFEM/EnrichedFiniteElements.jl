using Revise
using LinearAlgebra
using Test


using EnrichedFiniteElements



const Transformations = EnrichedFiniteElements.transformationFunctions
const integrator = EnrichedFiniteElements.Operators

const mesh_create = EnrichedFiniteElements.MeshCreation
const wave_func = EnrichedFiniteElements.EnrichmentCreator
const matrix_comp = EnrichedFiniteElements.MatrixCreation

function setup_test_environment(; wave_x::Int64 = 0, wave_y::Int64 = 0)
    #! Set up for the matrix creation
    domain = ((0, 1), (0, 1))

    # Create the mesh
    mesh = mesh_create.rectangle_domain(domain)

    nodes = mesh.nodes
    connectivity = mesh.connectivity
    boundary_index = mesh.boundary_idx
    boundary_edges = mesh.boundary_edges

    ansatz_wave = wave_func.create_wavenumbers(wave_x, wave_y)
    test_ansatz = wave_func.create_wavenumbers(wave_x, wave_y)
    wavenumbers_ansatz, wavenumbers_test =
        wave_func.wavenumber_creation(ansatz_wave, test_ansatz, 2)

    idx_wave_ansatz = collect(1:size(wavenumbers_ansatz, 1))
    idx_wave_test = collect(1:size(wavenumbers_test, 1))

    all_pairs = wave_func.generate_pairs(idx_wave_ansatz, idx_wave_test)
    idx_connectivity = collect(1:size(connectivity, 1))
    wave_node_pairs = wave_func.combine_wavenumber_with_all_nodes(all_pairs, connectivity)


    return nodes,
    connectivity,
    boundary_index,
    boundary_edges,
    wavenumbers_ansatz,
    wavenumbers_test,
    idx_wave_ansatz,
    idx_wave_test,
    all_pairs,
    idx_connectivity,
    wave_node_pairs
end

@testset "mass_matrix_jump_NO_ENRICH" begin
    dt = 0.1
    t0 = 0.0
    nodes,
    connectivity,
    boundary_index,
    boundary_edges,
    wavenumbers_ansatz,
    wavenumbers_test,
    idx_wave_ansatz,
    idx_wave_test,
    all_pairs,
    idx_connectivity,
    wave_node_pairs = setup_test_environment(wave_x = 0, wave_y = 0)

    cell_sparse_zero_array = matrix_comp.compute_sparse_mass_matrix(
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

end

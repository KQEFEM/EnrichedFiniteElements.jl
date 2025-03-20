using Revise
using LinearAlgebra
using Test


using EnrichedFiniteElements



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
@testset "Matrix Creation" begin

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

        array = matrix_comp.compute_sparse_mass_matrix(
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

        matrix = [
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

        @test isapprox(norm(real(array[1] - matrix)), 0.00022378579972373314, atol = 1e-3) # This is with the spatial step, dx
        @test isapprox(norm(imag(array[1])), 0.0, atol = 1e-12) # Ensures that the standard FEM method doesn't produce imaginary numbers 

    end

    @testset "mass_jump_enriched" begin
        @test isequal(1, 1)
    end


end

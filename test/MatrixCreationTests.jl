using Revise
using LinearAlgebra
using Test
using DelimitedFiles
using SparseArrays


using EnrichedFiniteElements



const integrator = EnrichedFiniteElements.Operators

const mesh_create = EnrichedFiniteElements.MeshCreation
const wave_func = EnrichedFiniteElements.EnrichmentCreator
const matrix_comp = EnrichedFiniteElements.MatrixCreation
#! Set up for the matrix creation
const domain = ((0, 1), (0, 1))
const mesh = mesh_create.rectangle_domain(domain)

function load_matlab_matrix(filename)
    # Load the data from the text file
    data = readdlm(filename)

    # Extract row indices, column indices, and values
    rows = Int.(data[:, 1])
    cols = Int.(data[:, 2])
    vals = complex.(data[:, 3], data[:, 4])  # Combine real and imaginary parts



    return sparse(rows, cols, vals)
end

function setup_test_environment(;
    wave_x::Int64 = 0,
    wave_y::Int64 = 0,
    zero_frequencies::Bool = false,
    time_enrichment_only::Bool = false,
)
    """
    zero_frequencies : bool: sets the frequencies to be 0
    """


    nodes = mesh.nodes
    connectivity = mesh.connectivity
    boundary_index = mesh.boundary_idx
    boundary_edges = mesh.boundary_edges

    ansatz_wave = wave_func.create_wavenumbers(wave_x, wave_y)
    test_ansatz = wave_func.create_wavenumbers(wave_x, wave_y)
    if zero_frequencies == false
        wavenumbers_ansatz, wavenumbers_test = wave_func.wavenumber_creation(
            ansatz_wave,
            test_ansatz,
            Number_of_Frequencies_per_Wavenumber = 2,
            time_enrichment_only = time_enrichment_only,
        )
    else
        wavenumbers_ansatz, wavenumbers_test = wave_func.wavenumber_creation(
            ansatz_wave,
            test_ansatz,
            Number_of_Frequencies_per_Wavenumber = 1,
            time_enrichment_only = time_enrichment_only,
        )
    end

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
    @testset "Mass Matrix" begin
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

            cell, _ = matrix_comp.compute_sparse_matrix(
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
                convection_bool = false,
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

            array = matrix_comp.convert_sparse_cell_to_array(cell[1])

            @test isapprox(norm(real(array - matrix)), 0.00022378579972373314, atol = 1e-3) # This is with the spatial step, dx
            @test isapprox(norm(imag(array[1])), 0.0, atol = 1e-12) # Ensures that the standard FEM method doesn't produce imaginary numbers 

        end

        @testset "mass_jump_enriched_noFrequencies" begin
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
            wave_node_pairs =
                setup_test_environment(wave_x = 1, wave_y = 1, zero_frequencies = true)

            cell, _ = matrix_comp.compute_sparse_matrix(
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
                convection_bool = false,
            )


            array = matrix_comp.convert_sparse_cell_to_array(cell[1])
            exact_matrix = conj(
                load_matlab_matrix("test/testdata/MassMatrixEnriched_noFrequencies.txt"),
            )
            #! This conjudate is simply as the matlab code orders in a differnet way

            @test isapprox(norm(array - exact_matrix), 3.984715840345388e-7) # This is with the spatial step, dx
            @test isapprox(norm(real(array - exact_matrix)), 2.92067550635342e-7)
            @test isapprox(norm(imag(array - exact_matrix)), 2.7106484307055905e-7)

        end

        @testset "mass_jump_enriched" begin
            """ Space-time enriched test"""

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
            wave_node_pairs =
                setup_test_environment(wave_x = 1, wave_y = 1, zero_frequencies = false)

            cell, _ = matrix_comp.compute_sparse_matrix(
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
                convection_bool = false,
            )

            array = matrix_comp.convert_sparse_cell_to_array(cell[1])

            exact_matrix =
                conj(load_matlab_matrix("test/testdata/MassMatrixEnriched_enriched.txt")) #! This conjudate is simply as the matlab code orders in a differnet way
            @test isapprox(norm(array - exact_matrix), 7.536714627244855e-7) # This is with the spatial step, dx
            @test isapprox(norm(real(array - exact_matrix)), 5.478068933256133e-7)
            @test isapprox(norm(imag(array - exact_matrix)), 5.176178912578366e-7)

        end
    end

    @testset "Convection Matrix" begin

        @testset "pDtq FEM" begin
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
            _, convection_cell = matrix_comp.compute_sparse_matrix(
                all_pairs,
                nodes,
                wave_node_pairs,
                wavenumbers_ansatz,
                wavenumbers_test,
                integrator,
                dt,
                convection_bool = true,
            )
            pDtq_cell = convection_cell[1]
            pDtq_cell = permutedims(pDtq_cell, (2, 1)) #!This fixes the ordering of the enrichments for julia to match MATLAB

            array = matrix_comp.convert_sparse_cell_to_array(pDtq_cell)

            @test isequal(array, spzeros(size(array, 1), size(array, 2)))

        end
        @testset "pDtq EFEM Time only" begin
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
            wave_node_pairs =
                setup_test_environment(wave_x = 1, wave_y = 1, time_enrichment_only = true)

            _, convection_cell = matrix_comp.compute_sparse_matrix(
                all_pairs,
                nodes,
                wave_node_pairs,
                wavenumbers_ansatz,
                wavenumbers_test,
                integrator,
                dt,
                convection_bool = true,
            )
            pDtq_cell = convection_cell[1]
            pDtq_cell = permutedims(pDtq_cell, (2, 1)) #!This fixes the ordering of the enrichments for julia to match MATLAB

            array = matrix_comp.convert_sparse_cell_to_array(pDtq_cell)
            exact_matrix = (load_matlab_matrix("test/testdata/ConvectionDt_time.txt"))
            @test isapprox(norm((array - exact_matrix)), 2.3595901292942768e-8)

            @test isapprox(norm(imag(array - exact_matrix)), 2.3511831062368027e-8)

            @test isapprox(norm(real(array - exact_matrix)), 1.9900648032073808e-9)


        end
        @testset "pDtq EFEM Space-Time" begin
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
            wave_node_pairs =
                setup_test_environment(wave_x = 1, wave_y = 1, time_enrichment_only = false)

            _, convection_cell = matrix_comp.compute_sparse_matrix(
                all_pairs,
                nodes,
                wave_node_pairs,
                wavenumbers_ansatz,
                wavenumbers_test,
                integrator,
                dt,
                convection_bool = true,
            )

            pDtq_cell = convection_cell[1]
            pDtq_cell = permutedims(pDtq_cell, (2, 1)) #!This fixes the ordering of the enrichments for julia to match MATLAB

            array = matrix_comp.convert_sparse_cell_to_array(pDtq_cell)
            exact_matrix = (load_matlab_matrix("test/testdata/ConvectionDt_space_time.txt"))
            @test isapprox(norm((array - exact_matrix)), 8.945366854027598e-8)

            @test isapprox(norm(imag(array - exact_matrix)), 6.459786881101117e-8)

            @test isapprox(norm(real(array - exact_matrix)), 6.187951325268306e-8)

        end

        @testset "v_nabla_q FEM" begin
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
            wave_node_pairs =
                setup_test_environment(wave_x = 0, wave_y = 0, time_enrichment_only = false)

            _, convection_cell = matrix_comp.compute_sparse_matrix(
                all_pairs,
                nodes,
                wave_node_pairs,
                wavenumbers_ansatz,
                wavenumbers_test,
                integrator,
                dt,
                convection_bool = true,
            )

            vDxeta_cell = convection_cell[2]

            vDxeta_cell = permutedims(vDxeta_cell, (2, 1)) #!There is a missing transpose somewhere in the basis operations
            array = matrix_comp.convert_sparse_cell_to_array(vDxeta_cell)

            exact_matrix = (load_matlab_matrix("test/testdata/ConvectionDX_FEM.txt"))
            @test isapprox(norm((array - exact_matrix)), 1.3442648471406098e-13)

            @test isapprox(norm(real(array - exact_matrix)), 1.3442648471406098e-13)

            @test isequal(norm(imag(array - exact_matrix)), 0.0)


        end
        @testset "v_nabla_q Time Enriched" begin
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
            wave_node_pairs =
                setup_test_environment(wave_x = 1, wave_y = 1, time_enrichment_only = true)

            _, convection_cell = matrix_comp.compute_sparse_matrix(
                all_pairs,
                nodes,
                wave_node_pairs,
                wavenumbers_ansatz,
                wavenumbers_test,
                integrator,
                dt,
                convection_bool = true,
            )

            vDxeta_cell = convection_cell[2]

            vDxeta_cell = permutedims(vDxeta_cell, (2, 1)) #!There is a missing transpose somewhere in the basis operations
            array = matrix_comp.convert_sparse_cell_to_array(vDxeta_cell)

            exact_matrix = (load_matlab_matrix("test/testdata/ConvectionDX_time.txt"))

            @test isapprox(norm((array - exact_matrix)), 6.715278246898072e-13)

            @test isapprox(norm(real(array - exact_matrix)), 6.694169592754642e-13)

            @test isapprox(norm(imag(array - exact_matrix)), 5.320295073598745e-14)
            @test isapprox(norm(diag(array - exact_matrix)), 1.285253769296312e-13)



        end

        @testset "v_nabla_q Space-Time Enriched" begin
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
            wave_node_pairs =
                setup_test_environment(wave_x = 1, wave_y = 1, time_enrichment_only = false)

            _, convection_cell = matrix_comp.compute_sparse_matrix(
                all_pairs,
                nodes,
                wave_node_pairs,
                wavenumbers_ansatz,
                wavenumbers_test,
                integrator,
                dt,
                convection_bool = true,
            )

            vDxeta_cell = convection_cell[2]

            vDxeta_cell = permutedims(vDxeta_cell, (2, 1)) #!There is a missing transpose somewhere in the basis operations
            array = matrix_comp.convert_sparse_cell_to_array(vDxeta_cell)

            exact_matrix = (load_matlab_matrix("test/testdata/ConvectionDX_space_time.txt"))

            @test isapprox(norm((array - exact_matrix)), 8.594550234390262e-8)

            @test isapprox(norm(real(array - exact_matrix)), 6.211649103163213e-8)

            @test isapprox(norm(imag(array - exact_matrix)), 5.939840835462623e-8)
            @test isapprox(norm(diag(array - exact_matrix)), 1.154607276802081e-8)



        end
    end



end

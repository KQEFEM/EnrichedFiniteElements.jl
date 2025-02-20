using Revise
using EnrichedFiniteElements

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

ansatz_wave = wave_func.create_wavenumbers(1, 1)
test_ansatz = wave_func.create_wavenumbers(1, 1)
wavenumbers_ansatz, wavenumbers_test =
    wave_func.wavenumber_creation(ansatz_wave, test_ansatz, 2)

idx_wave_ansatz = collect(1:size(wavenumbers_ansatz, 1))
idx_wave_test = collect(1:size(wavenumbers_test, 1))

all_pairs = wave_func.generate_pairs(idx_wave_ansatz, idx_wave_test)
println(all_pairs)
idx_connectivity = collect(1:size(connectivity, 1))
result = wave_func.combine_wavenumber_with_all_nodes(all_pairs, connectivity)

#! Matrix Creation 
using Revise
using SparseArrays
const Transformations = EnrichedFiniteElements.transformationFunctions
const integrator = EnrichedFiniteElements.Operators
cell_sparse_zero_array = Array{Any}(undef, 3, 3) # Or Array{SparseMatrixCSC{Float64, Int64}}(undef, 3, 3) for type specificity
n = size(nodes, 1)
for i = 1:9
    cell_sparse_zero_array[i] = spzeros(n, n) # Create nxn sparse matrix of zeros
end

v_nabla_q = 0
for ii in result
    # ii = result[1]
    println(ii[1][1][1])
    wave_ansatz_loc = wavenumbers_ansatz[ii[1][1][1], :]
    wave_test_loc = wavenumbers_test[ii[1][1][2], :]
    println(wave_ansatz_loc)
    println(wave_test_loc)
    triangle_connectivity = ii[2]
    println(triangle_connectivity)
    triangle_nodes = nodes[ii[2], :]
    println(triangle_nodes)

    triangle_nodes, triangle_connectivity =
        Transformations.correct_triangle_orientation!(triangle_nodes, triangle_connectivity)

    kx_kkx = wave_ansatz_loc[1] - wave_test_loc[1]
    ky_kky = wave_ansatz_loc[2] - wave_test_loc[2]
    A =
        kx_kkx * (triangle_nodes[2, 1] - triangle_nodes[1, 1]) +
        ky_kky * (triangle_nodes[2, 2] - triangle_nodes[1, 2])
    B =
        kx_kkx * (triangle_nodes[3, 1] - triangle_nodes[1, 1]) +
        ky_kky * (triangle_nodes[3, 2] - triangle_nodes[1, 2])
    C = kx_kkx * triangle_nodes[1, 1] + ky_kky * triangle_nodes[1, 2]

    tri_area, ddx, ddy =
        Transformations.Gradients_Larson(triangle_nodes[:, 1], triangle_nodes[:, 2])

    println(ddx)

    grads_grads_dx = ddx * ddx'
    grads_grads_dy = ddy * ddy'
    println(grads_grads_dx)
    println(triangle_nodes)


    #! Integration stage now 

    upper_bounds = [1.0, 1.0, 1.0]
    lower_bounds = [-1.0, -1.0, 0.0]
    omega = [0.0, 0.0]
    v_nabla_q =
        integrator.v_nabla_q(upper_bounds, lower_bounds, A, B, C, omega, grads_grads_dx, 0)
    println(v_nabla_q)

    #     v_nabla_q(
    #     upper_bounds::Vector{Float64},
    #     lower_bounds::Vector{Float64},
    #     A_val::Float64,
    #     B_val::Float64,
    #     C_val::Float64,
    #     omega::Vector{Float64},
    #     grads_matrix::Matrix{Float64},
    #     K,
    # )
    break
end

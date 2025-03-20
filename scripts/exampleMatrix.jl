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

wave_x = 0
wave_y = 0
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

#! Matrix Creation 
using Revise
using LinearAlgebra  # Import the LinearAlgebra module
using SparseArrays
const Transformations = EnrichedFiniteElements.TransformationFunctions
const integrator = EnrichedFiniteElements.Operators
cell_sparse_zero_array =
    Array{SparseMatrixCSC{ComplexF64,Int64}}(undef, size(all_pairs, 1), size(all_pairs, 1)) # Explicitly specify ComplexF64
n = size(nodes, 1)
for i = 1:prod(size(cell_sparse_zero_array))
    cell_sparse_zero_array[i] = spzeros(n, n) # Create nxn sparse matrix of zeros
end

v_nabla_q = 0
for ii in result
    # ii = result[1]

    wave_ansatz_loc = wavenumbers_ansatz[ii[1][1][1], :]
    wave_test_loc = wavenumbers_test[ii[1][1][2], :]


    triangle_connectivity = ii[2]

    triangle_nodes = nodes[ii[2], :]

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

    # println(ddx)

    grads_grads_dx = ddx * ddx'
    # println(grads_grads_dx)
    # println(size(grads_grads_dx))
    grads_grads_dy = ddy * ddy'
    # println(grads_grads_dx)
    # println(triangle_nodes)


    #! Integration stage now 

    upper_bounds = [1.0, 1.0, 1.0]
    lower_bounds = [0.0, -1.0, -1.0]
    omega = [0.0, 0.0]
    v_nabla_q_loc, _ = integrator.v_nabla_q(
        upper_bounds,
        lower_bounds,
        A,
        B,
        C,
        omega,
        grads_grads_dx,
        0,
        tri_area,
    )
    # println(v_nabla_q_loc)
    # println(cell_sparse_zero_array[1,2][triangle_connectivity,triangle_connectivity])
    cell_sparse_zero_array[1, 1][triangle_connectivity, triangle_connectivity] =
        cell_sparse_zero_array[1, 1][triangle_connectivity, triangle_connectivity] +
        v_nabla_q_loc

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
    # break
end
println(cell_sparse_zero_array[1, 1])
real(cell_sparse_zero_array[1])
real(diag(cell_sparse_zero_array[1]))
# imag(cell_sparse_zero_array[1])
cell_sparse_zero_array[1, 1][boundary_index, boundary_index] .= 0;
real(cell_sparse_zero_array[1])
real(cell_sparse_zero_array[1]) / 0.4375^2


############ MASS MATRIX ####################
#! This is being used to debug the above 

using Revise
using LinearAlgebra  # Import the LinearAlgebra module
using SparseArrays
const Transformations = EnrichedFiniteElements.TransformationFunctions
const integrator = EnrichedFiniteElements.Operators
cell_sparse_zero_array =
    Array{SparseMatrixCSC{ComplexF64,Int64}}(undef, size(all_pairs, 1), size(all_pairs, 1)) # Explicitly specify ComplexF64
n = size(nodes, 1)
for i = 1:prod(size(cell_sparse_zero_array))
    cell_sparse_zero_array[i] = spzeros(n, n) # Create nxn sparse matrix of zeros
end

mass_mat_loc = 0
for ii in result
    #? This is correct
    # ii = result[1]

    wave_ansatz_loc = wavenumbers_ansatz[ii[1][1][1], :]
    wave_test_loc = wavenumbers_test[ii[1][1][2], :]


    triangle_connectivity = ii[2]

    triangle_nodes = nodes[ii[2], :]

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

    # println(ddx)

    grads_grads_dx = ddx * ddx'
    # println(grads_grads_dx)
    # println(size(grads_grads_dx))
    grads_grads_dy = ddy * ddy'
    # println(grads_grads_dx)
    # println(triangle_nodes)


    #! Integration stage now (this is incorrect or there is a wrong sign)

    upper_bounds = [1.0, 1.0, 1.0]
    lower_bounds = [-1.0, -1.0, 0.0]
    omega = [0.0, 0.0]
    dt = 0.1
    t_jump = 0.0
    t0 = 0.0

    mass_loc, _ = integrator.mass_jump(
        upper_bounds,
        lower_bounds,
        A,
        B,
        C,
        omega,
        t_jump,
        dt,
        t0,
        tri_area,
    )
    # println(v_nabla_q_loc)
    # println(cell_sparse_zero_array[1,2][triangle_connectivity,triangle_connectivity])
    cell_sparse_zero_array[1, 1][triangle_connectivity, triangle_connectivity] =
        cell_sparse_zero_array[1, 1][triangle_connectivity, triangle_connectivity] +
        mass_loc

    # upper_bounds::Vector{Float64},
    # lower_bounds::Vector{Float64},
    # A_val::Float64,
    # B_val::Float64,
    # C_val::Float64,
    # omega::Vector{Float64},
    # dt::Real,
    # t0::Real,
    # triangle_area::Float64

end
println(cell_sparse_zero_array[1, 1])
real(cell_sparse_zero_array[1])
real(diag(cell_sparse_zero_array[1]))
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

real(cell_sparse_zero_array[1] - matrix)
diag(real(cell_sparse_zero_array[1] - matrix))
norm(diag(real(cell_sparse_zero_array[1] - matrix)))


############ Enriched MASS MATRIX ####################
#! This is being used to debug the above 

using Revise
using EnrichedFiniteElements

using LinearAlgebra  # Import the LinearAlgebra module
using SparseArrays
const Transformations = EnrichedFiniteElements.TransformationFunctions
const integrator = EnrichedFiniteElements.Operators

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
wavenumbers_ansatz, wavenumbers_test =
    wave_func.wavenumber_creation(ansatz_wave, test_ansatz, 2)

idx_wave_ansatz = collect(1:size(wavenumbers_ansatz, 1))
idx_wave_test = collect(1:size(wavenumbers_test, 1))

all_pairs = wave_func.generate_pairs(idx_wave_ansatz, idx_wave_test)
println(all_pairs)
idx_connectivity = collect(1:size(connectivity, 1))
result = wave_func.combine_wavenumber_with_all_nodes(all_pairs, connectivity)

cell_sparse_zero_array =
    Array{SparseMatrixCSC{ComplexF64,Int64}}(undef, size(all_pairs, 1), size(all_pairs, 1)) # Explicitly specify ComplexF64
n = size(nodes, 1)
for i = 1:prod(size(cell_sparse_zero_array))
    cell_sparse_zero_array[i] = spzeros(n, n) # Create nxn sparse matrix of zeros
end

mass_mat_loc = 0
for ii in result
    #? This is correct
    # ii = result[1]


    wave_ansatz_loc = wavenumbers_ansatz[ii[1][1][1], :]
    wave_test_loc = wavenumbers_test[ii[1][1][2], :]
    cell_idx = ii[1][1] # this grabs the tuple ( - , - )

    triangle_connectivity = ii[2]

    triangle_nodes = nodes[ii[2], :]

    triangle_nodes, triangle_connectivity =
        Transformations.correct_triangle_orientation!(triangle_nodes, triangle_connectivity)

    kx_kkx = wave_ansatz_loc[1] - wave_test_loc[1]
    ky_kky = wave_ansatz_loc[2] - wave_test_loc[2]
    omega = [wave_ansatz_loc[3], wave_test_loc[3]]
    A =
        kx_kkx * (triangle_nodes[2, 1] - triangle_nodes[1, 1]) +
        ky_kky * (triangle_nodes[2, 2] - triangle_nodes[1, 2])
    B =
        kx_kkx * (triangle_nodes[3, 1] - triangle_nodes[1, 1]) +
        ky_kky * (triangle_nodes[3, 2] - triangle_nodes[1, 2])
    C = kx_kkx * triangle_nodes[1, 1] + ky_kky * triangle_nodes[1, 2]

    tri_area, ddx, ddy =
        Transformations.Gradients_Larson(triangle_nodes[:, 1], triangle_nodes[:, 2])

    # println(ddx)

    grads_grads_dx = ddx * ddx'
    # println(grads_grads_dx)
    # println(size(grads_grads_dx))
    grads_grads_dy = ddy * ddy'
    # println(grads_grads_dx)
    # println(triangle_nodes)


    #! Integration stage now (this is incorrect or there is a wrong sign)

    upper_bounds = [1.0, 1.0, 1.0]
    lower_bounds = [-1.0, -1.0, 0.0]
    # omega = [0.0, 0.0]
    dt = 0.1
    t_jump = 0.0
    t0 = 0.0

    mass_loc, _ = integrator.mass_jump(
        upper_bounds,
        lower_bounds,
        A,
        B,
        C,
        omega,
        t_jump,
        dt,
        t0,
        tri_area,
    )
    # println(v_nabla_q_loc)
    # println(cell_sparse_zero_array[1,2][triangle_connectivity,triangle_connectivity])
    cell_sparse_zero_array[cell_idx[1], cell_idx[2]][
        triangle_connectivity,
        triangle_connectivity,
    ] .+= mass_loc

    # upper_bounds::Vector{Float64},
    # lower_bounds::Vector{Float64},
    # A_val::Float64,
    # B_val::Float64,
    # C_val::Float64,
    # omega::Vector{Float64},
    # dt::Real,
    # t0::Real,
    # triangle_area::Float64

end
println(cell_sparse_zero_array[1, 1])
real(cell_sparse_zero_array[1])
real(diag(cell_sparse_zero_array[1]))

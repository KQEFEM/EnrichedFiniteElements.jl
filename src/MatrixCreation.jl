module MatrixCreation
using SparseArrays

using HCubature
import ..BasisFunctions as basis
import ..Operators as op
# import ..TransformationFunctions as transformations
using ..TransformationFunctions

export all
const transformations = TransformationFunctions
const integrator = op

function sparse_matrix_creation(all_pairs, nodes)
    """
    Creates the sparse array needed fro the enriched matrix computations 
    """
    cell_sparse_zero_array = Array{SparseMatrixCSC{ComplexF64,Int64}}(
        undef,
        all_pairs[end][1],
        all_pairs[end][2],
    )
    n = size(nodes, 1)

    for i = 1:prod(size(cell_sparse_zero_array))
        cell_sparse_zero_array[i] = spzeros(n, n)
    end

    return cell_sparse_zero_array
end






function nodal_transformations(ii, nodes)
    """ Extracts the mesh values """
    triangle_connectivity = ii[2]
    triangle_nodes = nodes[ii[2], :]

    triangle_nodes, triangle_connectivity =
        transformations.correct_triangle_orientation!(triangle_nodes, triangle_connectivity)

    return triangle_nodes, triangle_connectivity
end

function enrichment_transformations(
    ii,
    wavenumbers_ansatz,
    wavenumbers_test,
    triangle_nodes,
)
    """ 
    This extracts the relavant variables for the enrichments and TransformationFunctions
    """
    wave_ansatz_loc = wavenumbers_ansatz[ii[1][1][1], :]
    wave_test_loc = wavenumbers_test[ii[1][1][2], :]
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
    return wave_ansatz_loc, wave_test_loc, kx_kkx, ky_kky, omega, A, B, C
end

function compute_sparse_matrix(
    all_pairs,
    nodes,
    connectivity_matrix,
    wavenumbers_ansatz,
    wavenumbers_test,
    integrator,
    dt; #? This allows me to use keyword arguments 
    t_jump = 0.0,
    t0 = 0.0,
    mass_bool::Bool = false,
    convection_bool::Bool = false,
)
    """
    This will create the matrices based on the booleon's given. 
    Note; All matrices will have the same size. 
    Note: This code currently needs to be ran separately for the mass jumps

    test_bool: will create the mass matrix 
    connectivity_matrix: joint connectivity with wavenumbers and nodes

    """
    indexing_array = Array{SparseMatrixCSC{ComplexF64,Int64}}(
        undef,
        all_pairs[end][1],
        all_pairs[end][2],
    ) #! This is simply for the for loop
    mass_sparse_array = nothing # Initialize outside the if block
    pDtq_array = nothing      # Initialize outside the if block
    vDxeta_array = nothing
    if mass_bool == true
        mass_sparse_array = sparse_matrix_creation(all_pairs, nodes)
    end
    if convection_bool == true
        pDtq_array = sparse_matrix_creation(all_pairs, nodes)
        #! There will be tohers to add also
        vDxeta_array = sparse_matrix_creation(all_pairs, nodes)

    end


    @views for (idx, ii) in enumerate(connectivity_matrix)
        #! This can be parallelised but it needs splitting into unique columns wrt to wave-pair 

        triangle_nodes, triangle_connectivity = nodal_transformations(ii, nodes)
        _, wave_test_loc , _, _, omega, A, B, C = enrichment_transformations(
            ii,
            wavenumbers_ansatz,
            wavenumbers_test,
            triangle_nodes,
        )

        tri_area, ddx, ddy =
            transformations.Gradients_Larson(triangle_nodes[:, 1], triangle_nodes[:, 2])
        grads_grads_dx = ddx .^ 2
        grads_grads_dy = ddy .^ 2


        cell_idx = LinearIndices(indexing_array)[ii[1][1][1], ii[1][1][2]] # this grabs the tuple ( - , - )

        if mass_bool == true
            upper_bounds = [1.0, 1.0, 1.0] # time integral is a dummy variable here
            lower_bounds = [-1.0, -1.0, 0.0] # time integral is a dummy variable here


            mass_sparse_array = create_components_mass_matrix(
                mass_sparse_array,
                cell_idx,
                triangle_connectivity,
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
        end
        if convection_bool == true
            # println(size(ddx))
            # println(ddx)
            pDtq_array, vDxeta_array =  create_components_convection_matrix(
            pDtq_array,
            vDxeta_array,
            cell_idx,
            triangle_connectivity,
            A,
            B,
            C,
            omega,
            dt,
            t0,
            ddx,
            wave_test_loc,
            tri_area,
        )
        pDtq_array = permutedims(pDtq_array, (2,1))
        vDxeta_array = permutedims(vDxeta_array, (2,1)) #!There is a missing transpose somewhere in the basis operations
        end
    end
    # mass_sparse_array = reshape(mass_sparse_array,sqrt(idx),:)
    return mass_sparse_array, pDtq_array, vDxeta_array
end

function convert_sparse_cell_to_array(
    sparse_cell_array::Matrix{SparseMatrixCSC{ComplexF64,Int64}},
)
    num_rows = size(sparse_cell_array, 1)
    num_cols = size(sparse_cell_array, 2)
    return reduce(
        vcat,
        [reduce(hcat, [sparse_cell_array[i, j] for j = 1:num_cols]) for i = 1:num_rows],
    )


end

function create_components_mass_matrix(
    mass_sparse_array::Array{SparseMatrixCSC{ComplexF64,Int64}},
    cell_idx::Int,
    triangle_connectivity::SubArray{
        Int64,
        1,
        Matrix{Int64},
        Tuple{Int64,Base.Slice{Base.OneTo{Int64}}},
        true,
    },
    upper_bounds::Vector{Float64},
    lower_bounds::Vector{Float64},
    A::Float64,
    B::Float64,
    C::Float64,
    omega::Vector{Float64},
    t_jump::Float64,
    dt::Float64,
    t0::Float64,
    tri_area::Float64,
)::Array{SparseMatrixCSC{ComplexF64,Int64}}
    """ 
    This creates the components for the mass matrix. It does not assemble mass matrix. 
    """
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
    mass_sparse_array[cell_idx][triangle_connectivity, triangle_connectivity] .+= mass_loc
    return mass_sparse_array
end

function create_components_convection_matrix(
    pDtq_array::Array{SparseMatrixCSC{ComplexF64,Int64}},
    vDxeta_array::Array{SparseMatrixCSC{ComplexF64,Int64}},
    cell_idx::Int,
    triangle_connectivity::SubArray{
        Int64,
        1,
        Matrix{Int64},
        Tuple{Int64,Base.Slice{Base.OneTo{Int64}}},
        true,
    },
    A::Float64,
    B::Float64,
    C::Float64,
    omega::Vector{Float64},
    dt::Float64,
    t0::Float64,
    ddx::Matrix{Float64},
    test_wavenumber::Vector{Float64},
    tri_area::Float64,
)#::(Array{SparseMatrixCSC{ComplexF64,Int64}},Array{SparseMatrixCSC{ComplexF64,Int64}})
"""
test_wavenumber
"""
            upper_bounds = [1.0, 1.0, dt]
            lower_bounds = [-1.0, -1.0, 0.0]


            pdtq_loc, _ =
                integrator.pDtq(upper_bounds, lower_bounds, A, B, C, omega, t0, tri_area)

            pDtq_array[cell_idx][triangle_connectivity, triangle_connectivity] .+= pdtq_loc
            vDxeta_loc, _ = integrator.v_nabla_q(upper_bounds, lower_bounds, A, B, C, omega, ddx, test_wavenumber[1], t0, tri_area)
            vDxeta_array[cell_idx][triangle_connectivity, triangle_connectivity] .+= vDxeta_loc

            return pDtq_array, vDxeta_array
end
end

#! # For check the matlab codes
# dt =0.1
# nodes = [ 0.0       0.0
# 1.0       0.0
# 1.0       1.0
# 0.0       1.0
# 0.5       0.0
# 1.0       0.5
# 0.5       1.0
# 0.0       0.5
# 0.29375   0.70625
# 0.375     0.375
# 0.647917  0.64375
# 0.71875   0.28125]
# nodes = [nodes, [1:size(nodes,1)].']
# connectivity = [  6   3  11
# 8   1  10
# 1   5  10
# 3   7  11
# 5   2  12
# 2   6  12
# 7   4   9
# 4   8   9
# 6  11  12
# 9   8  10
# 9  10  11
# 10   5  12
# 11  10  12
# 7   9  11]
# Connectivities = {Connectivity_Wavenumbers, connectivity,...
# connectivity};
#  Original_nodes_index = nodes, Periodic_Nodes = nodes

# matrix_connectivity_matrix = [
#     0.0312 0 0 0 0.0078 0 0 0.0078 0 0.0156 0 0
#     0 0.0234 0 0 0.0059 0.0059 0 0 0 0 0 0.0117
#     0 0 0.0295 0 0 0.0073 0.0074 0 0 0 0.0148 0
#     0 0 0 0.0245 0 0 0.0061 0.0061 0.0122 0 0 0
#     0.0078 0.0059 0 0 0.0371 0 0 0 0 0.0127 0 0.0107
#     0 0.0059 0.0073 0 0 0.0362 0 0 0 0 0.0122 0.0108
#     0 0 0.0074 0.0061 0 0 0.0368 0 0.0110 0 0.0123 0
#     0.0078 0 0 0.0061 0 0 0 0.0374 0.0109 0.0126 0 0
#     0 0 0 0.0122 0 0 0.0110 0.0109 0.0531 0.0094 0.0095 0
#     0.0156 0 0 0 0.0127 0 0 0.0126 0.0094 0.0697 0.0096 0.0098
#     0 0 0.0148 0 0 0.0122 0.0123 0 0.0095 0.0096 0.0682 0.0098
#     0 0.0117 0 0 0.0107 0.0108 0 0 0 0.0098 0.0098 0.0528
# ]

# # x_enrichments= 1
# # y_enrichments = 1
# # wavenumbers = reshape(collect(
# #         Iterators.product(-x_enrichments:x_enrichments, -y_enrichments:y_enrichments),
# #     ),:,1)


# # flattened_matrix = hcat(first.(wavenumbers), last.(wavenumbers))

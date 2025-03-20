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
        size(all_pairs, 1),
        size(all_pairs, 1),
    )
    n = size(nodes, 1)

    for i = 1:prod(size(cell_sparse_zero_array))
        cell_sparse_zero_array[i] = spzeros(n, n)
    end

    return cell_sparse_zero_array
end

# function mass_matrix(
#     nodes::Matrix{Float64},
#     connectivity::Matrix{Float64},
#     time,
#     wavenumbers::Matrix{Float64} = [0, 0, 0],
# )


# end






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
    return kx_kkx, ky_kky, omega, A, B, C
end

function compute_sparse_mass_matrix(
    all_pairs,
    nodes,
    result,
    wavenumbers_ansatz,
    wavenumbers_test,
    integrator,
    dt,
    t_jump = 0.0,
    t0 = 0.0,
)
    cell_sparse_zero_array = sparse_matrix_creation(all_pairs, nodes)

    for ii in result
        cell_idx = ii[1][1] # this grabs the tuple ( - , - )

        triangle_nodes, triangle_connectivity = nodal_transformations(ii, nodes)

        kx_kkx, ky_kky, omega, A, B, C = enrichment_transformations(
            ii,
            wavenumbers_ansatz,
            wavenumbers_test,
            triangle_nodes,
        )

        tri_area, ddx, ddy = transformations.Gradients_Larson(triangle_nodes[:, 1], triangle_nodes[:, 2])
        grads_grads_dx = ddx .^ 2
        grads_grads_dy = ddy .^ 2

        upper_bounds = [1.0, 1.0, 1.0]
        lower_bounds = [-1.0, -1.0, 0.0]

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

        cell_sparse_zero_array[cell_idx[1], cell_idx[2]][triangle_connectivity, triangle_connectivity] .+= mass_loc
    end

    return cell_sparse_zero_array
end

end


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

# matrix_result = [
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

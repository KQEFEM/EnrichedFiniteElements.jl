module MatrixCreation

using HCubature
import ..BasisFunctions as basis
import ..Operators as op
import ..transformationFunctions as Transformations

# const Transformations = .transformationFunctions
const integrator = op

function spares_matrix_creation(all_pairs,nodes)
    """
    Creates the sparse array needed fro the enriched matrix computations 
    """
    cell_sparse_zero_array = Array{SparseMatrixCSC{ComplexF64,Int64}}(undef, size(all_pairs, 1), size(all_pairs, 1))
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




end

function nodal_transformations(ii)
    """ Extracts the mesh values """
    triangle_connectivity = ii[2]
    triangle_nodes = nodes[ii[2], :]
    
    triangle_nodes, triangle_connectivity = 
        Transformations.correct_triangle_orientation!(triangle_nodes, triangle_connectivity)
    
    return triangle_nodes, triangle_connectivity
end

function enrichment_transformations(ii)
    """ 
    This extracts the relavant variables for the enrichments and transformationFunctions
    """
    wave_ansatz_loc = wavenumbers_ansatz[ii[1][1][1], :]
    wave_test_loc = wavenumbers_test[ii[1][1][2], :]
    kx_kkx = wave_ansatz_loc[1] - wave_test_loc[1]
        ky_kky = wave_ansatz_loc[2] - wave_test_loc[2]
        omega = [wave_ansatz_loc[3],wave_test_loc[3]]

        
        A = kx_kkx * (triangle_nodes[2, 1] - triangle_nodes[1, 1]) +
            ky_kky * (triangle_nodes[2, 2] - triangle_nodes[1, 2])
        B = kx_kkx * (triangle_nodes[3, 1] - triangle_nodes[1, 1]) +
            ky_kky * (triangle_nodes[3, 2] - triangle_nodes[1, 2])
        C = kx_kkx * triangle_nodes[1, 1] + ky_kky * triangle_nodes[1, 2]
        return kx_kkx, ky_kky, omega, A,B,C
end

function compute_sparse_mass_matrix(all_pairs, nodes, result, wavenumbers_ansatz, wavenumbers_test, integrator)
    cell_sparse_zero_array = matrix_creation(all_pairs,nodes);

    for ii in result
        cell_idx = ii[1][1] # this grabs the tuple ( - , - )

        triangle_nodes, triangle_connectivity = nodal_transformations(ii)
        
        kx_kkx, ky_kky, omega, A,B,C = enrichment_transformations(ii)
        
        tri_area, ddx, ddy = Transformations.Gradients_Larson(triangle_nodes[:, 1], triangle_nodes[:, 2])
        grads_grads_dx = ddx * ddx'
        grads_grads_dy = ddy * ddy'
        
        upper_bounds = [1.0, 1.0, 1.0]
        lower_bounds = [-1.0, -1.0, 0.0]

        dt = 0.1 #! This needs to be in the function input
        t_jump = 0.0
        t0 = 0.0

        mass_loc, _ = integrator.mass_jump(
            upper_bounds, lower_bounds, A, B, C, omega, t_jump, dt, t0, tri_area
        )
        
        cell_sparse_zero_array[cell_idx[1], cell_idx[2]][triangle_connectivity, triangle_connectivity] .+=
        mass_loc
    end
    
    return cell_sparse_zero_array
end


nodes = [ 0.0       0.0
1.0       0.0
1.0       1.0
0.0       1.0
0.5       0.0
1.0       0.5
0.5       1.0
0.0       0.5
0.29375   0.70625
0.375     0.375
0.647917  0.64375
0.71875   0.28125]
connectivity = [  6   3  11
8   1  10
1   5  10
3   7  11
5   2  12
2   6  12
7   4   9
4   8   9
6  11  12
9   8  10
9  10  11
10   5  12
11  10  12
7   9  11]


x_enrichments= 1
y_enrichments = 1
wavenumbers = reshape(collect(
        Iterators.product(-x_enrichments:x_enrichments, -y_enrichments:y_enrichments),
    ),:,1)


    flattened_matrix = hcat(first.(wavenumbers), last.(wavenumbers))

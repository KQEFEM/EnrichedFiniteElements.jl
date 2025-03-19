module MatrixCreation

using HCubature
import ..BasisFunctions as basis
import ..Operators as op

const Transformations = EnrichedFiniteElements.transformationFunctions
const integrator = EnrichedFiniteElements.Operators

function matrix_creation(all_pairs,nodes)
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

function mass_matrix(
    nodes::Matrix{Float64},
    connectivity::Matrix{Float64},
    time,
    wavenumbers::Matrix{Float64} = [0, 0, 0],
)


end




end

function compute_sparse_mass_matrix(all_pairs, nodes, result, wavenumbers_ansatz, wavenumbers_test, integrator)
    cell_sparse_zero_array = matrix_creation(all_pairs,nodes);

    for ii in result
        wave_ansatz_loc = wavenumbers_ansatz[ii[1][1][1], :]
        wave_test_loc = wavenumbers_test[ii[1][1][2], :]
        
        triangle_connectivity = ii[2]
        triangle_nodes = nodes[ii[2], :]
        
        triangle_nodes, triangle_connectivity = 
            Transformations.correct_triangle_orientation!(triangle_nodes, triangle_connectivity)
        
        kx_kkx = wave_ansatz_loc[1] - wave_test_loc[1]
        ky_kky = wave_ansatz_loc[2] - wave_test_loc[2]
        
        A = kx_kkx * (triangle_nodes[2, 1] - triangle_nodes[1, 1]) +
            ky_kky * (triangle_nodes[2, 2] - triangle_nodes[1, 2])
        B = kx_kkx * (triangle_nodes[3, 1] - triangle_nodes[1, 1]) +
            ky_kky * (triangle_nodes[3, 2] - triangle_nodes[1, 2])
        C = kx_kkx * triangle_nodes[1, 1] + ky_kky * triangle_nodes[1, 2]
        
        tri_area, ddx, ddy = Transformations.Gradients_Larson(triangle_nodes[:, 1], triangle_nodes[:, 2])
        grads_grads_dx = ddx * ddx'
        grads_grads_dy = ddy * ddy'
        
        upper_bounds = [1.0, 1.0, 1.0]
        lower_bounds = [-1.0, -1.0, 0.0]
        omega = [0.0, 0.0]
        dt = 0.1
        t_jump = 0.0
        t0 = 0.0

        mass_loc, _ = integrator.mass_jump(
            upper_bounds, lower_bounds, A, B, C, omega, t_jump, dt, t0, tri_area
        )
        
        cell_sparse_zero_array[1, 1][triangle_connectivity, triangle_connectivity] .+= mass_loc
    end
    
    return cell_sparse_zero_array
end

module LoadTerms

using SparseArrays

import ..BasisFunctions as basis
import ..Operators as integrator
import ..IntegrationWrappers as help
using ..TransformationFunctions

export all
const transformations = TransformationFunctions

function exact_pressure_0_velocity(kx::Real, ky::Real,  w::Real)
"""
This is based on the exact solution P49 in "Space-Time Enriched Finite Elements For
 Acoustic and Elastodynamic Problems"

 p = sin(k * xx - w*t), v=   k' / w * [sin(k * xx - w*t), -cos(k*x - w*t)]

 with RHS 

 f1 = ∂ₜp+∇v = (ky^2*cos(kx*x + ky*y - t*w))/w - w*cos(kx*x + ky*y - t*w) - (kx^2*sin(kx*x + ky*y - t*w))/w
 f_2 = ∂ₜv+∇p = [kx*cos(kx*x + ky*y - t*w) + kx*sin(kx*x + ky*y - t*w),0]

 returns: Functions of the above 

 Note: There are no transformations here 
"""

                                   
∂ₜp∇v(x::Real, y::Real, t::Real) = (ky^2*cos(kx*x + ky*y - t*w))/w - w*cos(kx*x + ky*y - t*w) - (kx^2*sin(kx*x + ky*y - t*w))/w
∂ₜv∇ₓp(x::Real, y::Real, t::Real) = kx*cos(kx*x + ky*y - t*w) + kx*sin(kx*x + ky*y - t*w)
∂ₜv∇ᵥp(x::Real, y::Real, t::Real) = 0.0

    # Return them as a tuple of functions
    return ∂ₜp∇v, ∂ₜv∇ₓp, ∂ₜv∇ᵥp
end

function load_term_integration( all_pairs,
    nodes,
    connectivity_matrix,
    wavenumbers_ansatz,
    wavenumbers_test;
    dt,
    t0,
    rhs_function, 
    args...)

    indexing_array = Array{SparseMatrixCSC{ComplexF64,Int64}}(
        undef,
        all_pairs[end][1],
        1,
    ) #! This is simply for the for loop
    f1 = sparse_matrix_creation(all_pairs, nodes)
    # Initialize outside the if block
    f2x = f1
    f2y = f1
    println(rhs_function)
    _,_, f2y_term = rhs_function
    f1_term = rhs_function[1]
    f2x_term = rhs_function[2]
    rhs_function[1](1,2,3)  # Assign the function to a variable
    rhs_function[2](1,2,3)  # Assign the function to a variable
    rhs_function[3](1,2,3)
    @views for (idx, ii) in enumerate(connectivity_matrix)

        triangle_nodes, triangle_connectivity = nodal_transformations(ii, nodes)
        _, wave_test_loc, _, _, omega, A, B, C = enrichment_transformations(
            ii,
            wavenumbers_ansatz,
            wavenumbers_test,
            triangle_nodes,
        )

        tri_area, gradients =
            transformations.Gradients_Larson(triangle_nodes[:, 1], triangle_nodes[:, 2])


        cell_idx = LinearIndices(indexing_array)[ii[1][1][1], ii[1][1][2]] # this grabs the tuple ( - , - )
        upper_bounds = [1.0, 1.0, dt]
        lower_bounds = [-1.0, -1.0, 0.0] #! this is wrong as it should be t_0,t_1
 # Integrating ∂ₜp
#  print(f1_term)
#  println(typeof(f1_term))
#  println(typeof(triangle_nodes))
 
 f1_integral, _ = integrator.load_term(upper_bounds, lower_bounds, A, B, C, omega, t0, tri_area,triangle_nodes, f1_term)
 
 # Integrating ∇ₓp
 f2x_integral , _= integrator.load_term(upper_bounds, lower_bounds, A, B, C, omega, t0, tri_area,triangle_nodes, f2x_term)
 
 # Integrating ∇ᵥp
 f2y_integral , _= integrator.load_term(upper_bounds, lower_bounds, A, B, C, omega, t0, tri_area,triangle_nodes, f2y_term)

#  println(f1_in print(dhsdjah)
# tegral)
#  println(f1)
#  println(size(f1_integral))
#  println(size(f1[cell_idx][triangle_connectivity]))
 f1[cell_idx][triangle_connectivity,1] .+= f1_integral
 f2x[cell_idx][triangle_connectivity,1] .+= f2x_integral
 f2y[cell_idx][triangle_connectivity,1] .+= f2y_integral

 
    end 
    return f1, f2x,f2y
end 


function sparse_matrix_creation(all_pairs, nodes)
    """
    Creates the sparse array needed fro the enriched matrix computations 
    """
    cell_sparse_zero_array = Array{SparseMatrixCSC{ComplexF64,Int64}}(
        undef,
        all_pairs[end][1],
        1,
    )
    n = size(nodes, 1)

    for i = 1:prod(size(cell_sparse_zero_array))
        cell_sparse_zero_array[i] = spzeros(n, 1)
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
    kx_kkx =  - wave_test_loc[1]
    ky_kky =  - wave_test_loc[2]
    omega =  - wave_test_loc[3]


    A =
        kx_kkx * (triangle_nodes[2, 1] - triangle_nodes[1, 1]) +
        ky_kky * (triangle_nodes[2, 2] - triangle_nodes[1, 2])
    B =
        kx_kkx * (triangle_nodes[3, 1] - triangle_nodes[1, 1]) +
        ky_kky * (triangle_nodes[3, 2] - triangle_nodes[1, 2])
    C = kx_kkx * triangle_nodes[1, 1] + ky_kky * triangle_nodes[1, 2]
    return wave_ansatz_loc, wave_test_loc, kx_kkx, ky_kky, omega, A, B, C
end


end 
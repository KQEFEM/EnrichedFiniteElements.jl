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

    returns:
    mass_cell: tuple containing the following
        mass_cell_t1: this is the upper time jump (t_1-t_0)
        mass_cell_t0: this is the lower time jump (t_0-t_0)
    
    convection_cell: tuple containing the following
        (pDtq_cell, vDxeta_cell,vDyeta_cell)

    """
    indexing_array = Array{SparseMatrixCSC{ComplexF64,Int64}}(
        undef,
        all_pairs[end][1],
        all_pairs[end][2],
    ) #! This is simply for the for loop
    mass_cell_t1 = nothing # Initialize outside the if block
    mass_cell_t0 = nothing
    pDtq_cell = nothing      # Initialize outside the if block
    vDxeta_cell = nothing
    vDyeta_cell = nothing
    if mass_bool == true
        mass_cell_t1 = sparse_matrix_creation(all_pairs, nodes)
        mass_cell_t0 = sparse_matrix_creation(all_pairs, nodes)
        mass_cell = nothing
    end
    if convection_bool == true
        pDtq_cell = sparse_matrix_creation(all_pairs, nodes)
        vDxeta_cell = sparse_matrix_creation(all_pairs, nodes)
        vDyeta_cell = sparse_matrix_creation(all_pairs, nodes)
        convection_cell = nothing

    end


    @views for (idx, ii) in enumerate(connectivity_matrix)
        #! This can be parallelised but it needs splitting into unique columns wrt to wave-pair 

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

        if mass_bool == true
            upper_bounds = [1.0, 1.0, 1.0] # time integral is a dummy variable here
            lower_bounds = [-1.0, -1.0, 0.0] # time integral is a dummy variable here


            mass_cell_t1 = create_components_mass_matrix(
                mass_cell_t1,
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
            mass_cell_t0 = create_components_mass_matrix(
                mass_cell_t0,
                cell_idx,
                triangle_connectivity,
                upper_bounds,
                lower_bounds,
                A,
                B,
                C,
                omega,
                t0,
                dt,
                t0,
                tri_area,
            )
        end
        if convection_bool == true
            # println(size(ddx))
            # println(ddx)
            pDtq_cell, vDxeta_cell, vDyeta_cell = create_components_convection_matrix(
                pDtq_cell,
                vDxeta_cell,
                vDyeta_cell,
                cell_idx,
                triangle_connectivity,
                A,
                B,
                C,
                omega,
                dt,
                t0,
                gradients,
                wave_test_loc,
                tri_area,
            )
        end
    end
    mass_cell = (mass_cell_t1, mass_cell_t0)
    convection_cell = (pDtq_cell, vDxeta_cell, vDyeta_cell)

    # mass_cell_t1 = reshape(mass_cell_t1,sqrt(idx),:)
    return mass_cell, convection_cell
end
function convert_sparse_cell_to_array(
    sparse_cell_array::Matrix{SparseMatrixCSC{ComplexF64,Int64}},
)
    """
    Converts a matrix of sparse matrices into a single dense matrix.
    
    This function takes a matrix where each entry is a sparse matrix and combines them
    into a single large matrix by horizontally concatenating the elements in each row
    and then vertically concatenating the resulting blocks.
    
    # Arguments
    - `sparse_cell_array::Matrix{SparseMatrixCSC{ComplexF64,Int64}}`: A matrix where
      each element is a sparse matrix.
    
    # Returns
    - A single `SparseMatrixCSC{ComplexF64,Int64}` combining all elements of `sparse_cell_array`.
    """
    num_rows = size(sparse_cell_array, 1)
    num_cols = size(sparse_cell_array, 2)
    return reduce(
        vcat,
        [reduce(hcat, [sparse_cell_array[i, j] for j = 1:num_cols]) for i = 1:num_rows],
    )
end

function create_components_mass_matrix(
    mass_cell::Array{SparseMatrixCSC{ComplexF64,Int64}},
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
    Computes and updates the components of the mass matrix.
    
    This function does not assemble the full mass matrix but updates the mass matrix
    components for a given cell index using integration over triangular elements.
    
    # Arguments
    - `mass_cell::Array{SparseMatrixCSC{ComplexF64,Int64}}`: Array of sparse mass matrices.
    - `cell_idx::Int`: Index of the current cell.
    - `triangle_connectivity::SubArray{Int64,1,...}`: Connectivity information for triangles.
    - `upper_bounds::Vector{Float64}`: Upper integration limits.
    - `lower_bounds::Vector{Float64}`: Lower integration limits.
    - `A::Float64, B::Float64, C::Float64`: Coefficients related to the mass matrix.
    - `omega::Vector{Float64}`: Frequency-related parameters.
    - `t_jump::Float64`: Jump discontinuity in time.
    - `dt::Float64`: Time step size.
    - `t0::Float64`: Initial time.
    - `tri_area::Float64`: Area of the triangle element.
    
    # Returns
    - `Array{SparseMatrixCSC{ComplexF64,Int64}}`: Updated mass matrix array.
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
    mass_cell[cell_idx][triangle_connectivity, triangle_connectivity] .+= mass_loc
    return mass_cell
end

function create_components_convection_matrix(
    pDtq_cell::Array{SparseMatrixCSC{ComplexF64,Int64}},
    vDxeta_cell::Array{SparseMatrixCSC{ComplexF64,Int64}},
    vDyeta_cell::Array{SparseMatrixCSC{ComplexF64,Int64}},
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
    gradients::Tuple{Matrix{Float64},Matrix{Float64}},
    test_wavenumber::Vector{Float64},
    tri_area::Float64,
)
    """
    Computes and updates components of the convection matrix.
    
    This function updates the convection-related matrices using integration over triangular elements,
    including contributions from time derivatives and spatial gradients.
    
    # Arguments
    - `pDtq_cell::Array{SparseMatrixCSC{ComplexF64,Int64}}`: Array of sparse matrices for time derivative terms.
    - `vDxeta_cell::Array{SparseMatrixCSC{ComplexF64,Int64}}`: Array of sparse matrices for x-gradient convection terms.
    - `vDyeta_cell::Array{SparseMatrixCSC{ComplexF64,Int64}}`: Array of sparse matrices for y-gradient convection terms.
    - `cell_idx::Int`: Index of the current cell.
    - `triangle_connectivity::SubArray{Int64,1,...}`: Connectivity information for triangles.
    - `A::Float64, B::Float64, C::Float64`: Coefficients used in integration.
    - `omega::Vector{Float64}`: Frequency-related parameters.
    - `dt::Float64`: Time step size.
    - `t0::Float64`: Initial time.
    - `gradients::Tuple{Matrix{Float64},Matrix{Float64}}`: Gradient matrices.
    - `test_wavenumber::Vector{Float64}`: Wavenumber parameters for testing.
    - `tri_area::Float64`: Area of the triangle element.
    
    # Returns
    - `(Array{SparseMatrixCSC{ComplexF64,Int64}}, Array{SparseMatrixCSC{ComplexF64,Int64}}, Array{SparseMatrixCSC{ComplexF64,Int64}})`: 
      Updated convection-related matrices for time, x-gradient, and y-gradient terms.
    """
    upper_bounds = [1.0, 1.0, dt]
    lower_bounds = [-1.0, -1.0, 0.0]
    ddx = gradients[1]
    ddy = gradients[2]

    pdtq_loc, _ = integrator.pDtq(upper_bounds, lower_bounds, A, B, C, omega, t0, tri_area)

    pDtq_cell[cell_idx][triangle_connectivity, triangle_connectivity] .+= pdtq_loc
    vDxeta_loc, _ = integrator.v_nabla_q(
        upper_bounds,
        lower_bounds,
        A,
        B,
        C,
        omega,
        ddx,
        test_wavenumber[1],
        t0,
        tri_area,
    )
    vDxeta_cell[cell_idx][triangle_connectivity, triangle_connectivity] .+= vDxeta_loc
    vDyeta_loc, _ = integrator.v_nabla_q(
        upper_bounds,
        lower_bounds,
        A,
        B,
        C,
        omega,
        ddy,
        test_wavenumber[1],
        t0,
        tri_area,
    )
    vDyeta_cell[cell_idx][triangle_connectivity, triangle_connectivity] .+= vDyeta_loc

    return pDtq_cell, vDxeta_cell, vDyeta_cell
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

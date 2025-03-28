module Formulation

using SparseArrays


function combine_matrices(;mass_cell = mass_cell, convection_cell = convection_cell)
    pDtq_array, vDxeta_array, vDyeta_array = convection_cell
    mass_array_t1, mass_array_t0 = mass_cell

    rows, cols = size(pDtq_array)
    # Create a sparse zero matrix of the same size
    sparse_zero_matrix = spzeros(rows, cols)

    zero_array = SparseMatrixCSC
    convection_matrix = [pDtq_array vDxeta_array vDyeta_array;
     vDxeta_array pDtq_array sparse_zero_matrix;
     vDyeta_array sparse_zero_matrix pDtq_array]

     mass_t1 = [mass_array_t1 sparse_zero_matrix sparse_zero_matrix;
     sparse_zero_matrix mass_array_t1 sparse_zero_matrix;
     sparse_zero_matrix sparse_zero_matrix mass_array_t1
     ]
     mass_t0 = [mass_array_t0 sparse_zero_matrix sparse_zero_matrix;
     sparse_zero_matrix mass_array_t0 sparse_zero_matrix;
     sparse_zero_matrix sparse_zero_matrix mass_array_t0
     ]

     return convection_matrix, mass_t1,mass_t0
end 




end 


module IntegrationWrappers
using ..BasisFunctions  # Or import BasisFunctions as basis
using ..TransformationFunctions

export all
const transformations = TransformationFunctions
export all

function time_enrichment_wrapper(v, f, omega, t0::Real)
    x = v[1]
    y = v[2]
    z = v[3]
    return f(x, y, z, omega, t0)
end

function space_enrichment_wrapper(v, f, A_val, B_val, C_val)
    x = v[1]
    y = v[2]
    z = v[3]
    return f(x, y, z, A_val, B_val, C_val)  # Call your original function
end
function hat_wrapper(v, f)
    x = v[1]
    y = v[2]
    z = v[3]
    return f(x, y, z)
end

function grads_wrapper(v, f, grad_matrix_input::AbstractMatrix{<:Real})
    """ Simply wraps the gradients matrix to allow for integration """

    x = v[1]
    y = v[2]

    z = v[3]
    grad_matrix = f(grad_matrix_input, x, y, z)  # Pass the matrix AND x, y, z

    return grad_matrix
end

function mass_jump_wrapper(v, f, w::Real, t_jump::Real, dt::Real, t0::Real, ww::Real)
    x = v[1]
    y = v[2]
    z = v[3]
    return f(x, y, z, w, t_jump, dt, t0, ww)
end

function grads_wrap(v)
    """ Simply wraps the gradients matrix to allow for integration """
    grads_wrapper(v, BasisFunctions.grads_matrix, grads_matrix)
end
function enrich_space_wrap(v)
    return space_enrichment_wrapper(v, BasisFunctions.enrich_space, A_val, B_val, C_val) # Or helper.space_enrichment_wrapper if aliased
end

function e_time_wrap(v)
    return time_enrichment_wrapper(v, BasisFunctions.enrichment_time, omega)
end

function hat_wrap(v)
    return hat_wrapper(v, BasisFunctions.p_mat)
end

function load_term_wrapper(v,f, triangle_nodes::Matrix{Float64})
    """ This should accept any number of arguments """
    x = v[1]
    y = v[2]
    z = v[3]
    # println(triangle_nodes)
    # println("Hello")
    X21,X31,Y21,Y31,X1,Y1 = transformations.arb_triangle_to_ref(triangle_nodes)
     # Apply transformation arb ▷ -> ref triangle ∟
     xx = X21 * x + X31 * y + X1
     yy = Y21 * x + Y31 * y + Y1
    return f(xx, yy, z)
end 

end

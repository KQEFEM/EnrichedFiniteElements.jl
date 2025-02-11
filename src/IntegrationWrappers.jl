module IntegrationWrappers
using ..BasisFunctions  # Or import BasisFunctions as basis

export all

function time_enrichment_wrapper(v, f, omega)
    x = v[1]
    y = v[2]
    z = v[3]
    return f(x, y, z, omega)
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

function mass_jump_wrapper(
    v,
    f,
    w::Vector{Float64},
    dt::Real,
    t0::Real,
    ww::Vector{Float64},
)
    x = v[1]
    y = v[2]
    z = v[3]
    return f(x, y, z, w, dt, t0, ww)
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


end

module Operators

using HCubature
import ..BasisFunctions as basis
import ..IntegrationWrappers as help

export all

"""
upper_bounds & lower_bounds: [dx,dy,dt]
omega: [w,ww]
K: [kx - kkx, ky - kky]
"""

function pDtq(
    upper_bounds::Vector{Float64},
    lower_bounds::Vector{Float64},
    A_val::Float64,
    B_val::Float64,
    C_val::Float64,
    omega::Vector{Float64},
)

    function integrand(v)
        x = v[1]
        y = v[2]
        z = v[3]
        space_enrichment = help.space_enrichment_wrapper(
            [(x + 1) / 2, (y + 1) / 2, z],
            basis.enrich_space,
            A_val,
            B_val,
            C_val,
        ) # Pass vector, call enrich_space
        time_enrichment = help.time_enrichment_wrapper(
            [(x + 1) / 2, (y + 1) / 2, z],
            basis.enrichment_time,
            diff(omega)[1],
        ) # Pass vector, call enrichment_time
        hat_ansatz = help.hat_wrapper([(x + 1) / 2, (y + 1) / 2, z], basis.phi) # Pass vector, call hat_function
        hat_test = help.hat_wrapper([(x + 1) / 2, (y + 1) / 2, z], basis.phi) # Pass vector, call hat_function

        return 1 / 4 *
               (1 - x) *
               space_enrichment *
               time_enrichment *
               hat_ansatz *
               (omega[2] * hat_test')  # Multiply the RESULTS
    end
    integral_result = hcubature(integrand, lower_bounds, upper_bounds) # Use integrand
    return integral_result
end

function v_nabla_q(
    upper_bounds::Vector{Float64},
    lower_bounds::Vector{Float64},
    A_val::Float64,
    B_val::Float64,
    C_val::Float64,
    omega::Vector{Float64},
    grads_matrix::Matrix{Float64},
    K,
)
    """ This is the term \\int_Omega v \\cdot \\nabla conj(q) and equivalent """
    function integrand(v)
        x = v[1]
        y = v[2]
        z = v[3]
        space_enrichment = help.space_enrichment_wrapper(
            [(x + 1) / 2, (y + 1) / 2, z],
            basis.enrich_space,
            A_val,
            B_val,
            C_val,
        ) # Pass vector, call enrich_space
        time_enrichment = help.time_enrichment_wrapper(
            [(x + 1) / 2, (y + 1) / 2, z],
            basis.enrichment_time,
            diff(omega)[1],
        ) # Pass vector, call enrichment_time
        hat_ansatz = help.hat_wrapper([(x + 1) / 2, (y + 1) / 2, z], basis.phi) # Pass vector, call hat_function
        hat_test = help.hat_wrapper([(x + 1) / 2, (y + 1) / 2, z], basis.phi) # Pass vector, call hat_function
        grads = help.grads_wrapper(
            [(x + 1) / 2, (y + 1) / 2, z],
            basis.grads_matrix,
            grads_matrix,
        )

        return 1 / 4 *
               (1 - x) *
               (hat_ansatz .* grads - 1im * hat_ansatz * hat_test' * K) *
               space_enrichment *
               time_enrichment
    end
    integral_result = hcubature(integrand, lower_bounds, upper_bounds) # Use integrand
    return integral_result
end

function mass_jump(
    upper_bounds::Vector{Float64},
    lower_bounds::Vector{Float64},
    A_val::Float64,
    B_val::Float64,
    C_val::Float64,
    omega::Vector{Float64},
    dt::Real,
    t0::Real,
)
    """ 2D integral, for example int_{d\\Omega_t} p^-[[ conj(q) ]]. This is the integral in 2D over the time boundary """
    function integrand(v)
        x = v[1]
        y = v[2]
        z = v[3]
        space_enrichment = help.space_enrichment_wrapper(
            [(x + 1) / 2, (y + 1) / 2, z],
            basis.enrich_space,
            A_val,
            B_val,
            C_val,
        ) # Pass vector, call enrich_space
        time_enrichment = help.mass_jump_wrapper(
            [(x + 1) / 2, (y + 1) / 2, z],
            basis.e_time_mass,
            omega[1],
            dt,
            t0,
            omega[2],
        )
        hat_ansatz = help.hat_wrapper([(x + 1) / 2, (y + 1) / 2, z], basis.phi) # Pass vector, call hat_function
        hat_test = help.hat_wrapper([(x + 1) / 2, (y + 1) / 2, z], basis.phi) # Pass vector, call hat_function

        return 1 / 4 *
               (1 - x) *
               (hat_ansatz * hat_test') *
               space_enrichment *
               time_enrichment
    end
    integral_result = hcubature(integrand, lower_bounds, upper_bounds) # Use integrand
    return integral_result
end

# function load_term()
#     function integrand(v)
#         x = v[1]
#         y = v[2]
#         z = v[3]
#         space_enrichment = help.space_enrichment_wrapper([(x + 1)/2, (y+1)/2, z], basis.enrich_space, A_val, B_val, C_val) # Pass vector, call enrich_space
#         time_enrichment = help.mass_jump_wrapper([(x + 1)/2, (y+1)/2, z],basis.e_time_mass, omega[1], dt, t0, omega[2])
#         hat_ansatz = help.hat_wrapper([(x + 1)/2, (y+1)/2, z], basis.phi) # Pass vector, call hat_function
#         hat_test = help.hat_wrapper([(x + 1)/2, (y+1)/2, z], basis.phi) # Pass vector, call hat_function

#         return (hat_ansatz *  hat_test') * space_enrichment * time_enrichment
#     end
#     integral_result = hcubature(integrand, lower_bounds, upper_bounds) # Use integrand
#     return integral_result
# end 
end # End of module

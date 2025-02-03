module Operators

using HCubature
import ..BasisFunctions as basis
import ..IntegrationWrappers as help

export all

function pDtq(upper_bounds::Vector{Float64}, lower_bounds::Vector{Float64}, A_val::Float64, B_val::Float64, C_val::Float64, omega::Vector{Float64}) 
    
    function integrand(v)
        x = v[1]
        y = v[2]
        z = v[3]
        space_enrichment = help.space_enrichment_wrapper([x, y, z], basis.enrich_space, A_val, B_val, C_val) # Pass vector, call enrich_space
        time_enrichment = help.time_enrichment_wrapper([x, y, z], basis.enrichment_time, diff(omega)[1]) # Pass vector, call enrichment_time
        hat_ansatz = help.hat_wrapper([x, y, z], basis.phi) # Pass vector, call hat_function
        hat_test = help.hat_wrapper([x, y, z], basis.phi) # Pass vector, call hat_function

        return space_enrichment * time_enrichment * hat_ansatz * (omega[2]  * hat_test')  # Multiply the RESULTS
    end
    integral_result = hcubature(integrand, lower_bounds, upper_bounds) # Use integrand
    return integral_result
end

    function v_nabla_q(upper_bounds::Vector{Float64}, lower_bounds::Vector{Float64}, A_val::Float64, B_val::Float64, C_val::Float64, omega::Vector{Float64}, grads_matrix::Matrix{Float64},K) 
    """ This is the term \\int_Omega v \\cdot \\nabla conj(q) and equivalent """
        function integrand(v)
            x = v[1]
            y = v[2]
            z = v[3]
            space_enrichment = help.space_enrichment_wrapper([x, y, z], basis.enrich_space, A_val, B_val, C_val) # Pass vector, call enrich_space
            time_enrichment = help.time_enrichment_wrapper([x, y, z], basis.enrichment_time, diff(omega)[1]) # Pass vector, call enrichment_time
            hat_ansatz = help.hat_wrapper([x, y, z], basis.phi) # Pass vector, call hat_function
            hat_test = help.hat_wrapper([x, y, z], basis.phi) # Pass vector, call hat_function
            grads = help.grads_wrapper([x,y,z], basis.grads_matrix,grads_matrix)

            return (hat_ansatz .* grads - 1im * hat_ansatz * hat_test' * K) * space_enrichment * time_enrichment
        end
        integral_result = hcubature(integrand, lower_bounds, upper_bounds) # Use integrand
        return integral_result

    end 
end # End of module
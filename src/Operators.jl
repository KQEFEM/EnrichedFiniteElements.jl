module Operators

using HCubature
import ..BasisFunctions as basis
import ..IntegrationWrappers as help

export all

function pDtq(upper_bounds::Vector{Float64}, lower_bounds::Vector{Float64}, A_val::Vector{Float64}, B_val::Vector{Float64}, C_val::Vector{Float64}, omega::Vector{Float64})    
    function integrand(v)
        x = v[1]
        y = v[2]
        z = v[3]
        space_result_ansatz = help.space_enrichment_wrapper([x, y, z], basis.e_space, A_val[1], B_val[1], C_val[1]) # Pass vector, call e_space
        time_result_ansatz = help.time_enrichment_wrapper([x, y, z], basis.e_time_ansatz, omega[1]) # Pass vector, call e_time_ansatz
        hat_result_ansatz = help.hat_wrapper([x, y, z], basis.phi) # Pass vector, call hat_function
        space_result_test = help.space_enrichment_wrapper([x, y, z], basis.e_space, A_val[2], B_val[2], C_val[2]) # Pass vector, call e_space
        time_result_test = help.time_enrichment_wrapper([x, y, z], basis.e_time_ansatz, omega[2]) # Pass vector, call e_time_ansatz
        hat_result_test = help.hat_wrapper([x, y, z], basis.phi) # Pass vector, call hat_function

        return space_result_ansatz * time_result_ansatz * hat_result_ansatz * conj(omega[2]space_result_test * time_result_test * hat_result_test')  # Multiply the RESULTS
    end
    integral_result = hcubature(integrand, lower_bounds, upper_bounds) # Use integrand
    return integral_result
end

end # End of module
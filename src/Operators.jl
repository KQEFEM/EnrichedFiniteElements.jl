module Operators

using HCubature
import ..BasisFunctions as basis
import ..IntegrationWrappers as help

export all

function pDtq(upper_bounds, lower_bounds, A_val, B_val, C_val, omega) # Add omega as argument
    function integrand(v)
        x = v[1]
        y = v[2]
        z = v[3]
        space_result = help.space_enrichment_wrapper([x, y, z], basis.e_space, A_val, B_val, C_val) # Pass vector, call e_space
        time_result = help.time_enrichment_wrapper([x, y, z], basis.e_time_ansatz, omega) # Pass vector, call e_time_ansatz
        hat_result = help.hat_wrapper([x, y, z], basis.phi) # Pass vector, call hat_function
        return space_result * time_result * hat_result  # Multiply the RESULTS
    end
    integral_result = hcubature(integrand, lower_bounds, upper_bounds) # Use integrand
    return integral_result
end

end # End of module
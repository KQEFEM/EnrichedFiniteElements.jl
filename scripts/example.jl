using Revise  # If you're using Revise (recommended)
using EnrichedFiniteElements

## Testing the operators 

EnrichedFiniteElements.BasisFunctions.phi(0.5, 0.5) # Test a function
println("Test successful!") # To confirm that the code ran

# Example of shorter alias:
const basis = EnrichedFiniteElements.BasisFunctions
basis.p_matrix(0.5, 0.5)
println("Alias Test successful!") # To confirm that the code ran

grad_matrix_input = [1.0 2.0 3.0; 4.0 5.0 6.0; 7.0 8.0 9.0] # Example matrix
grad_matrix = basis.grads_matrix(grad_matrix_input)


## Testing the integration 

upper_bounds = [1.0, 1.0, 1.0]
lower_bounds = [0.0, 0.0, 0.0]
A = 0.0
B = 0.0
C = 0.0
omega = [0.0,0] 
# A = 1.0 
# B = 0.2 
# C = 0.8
# omega = [3.0, - 3]
const operators = EnrichedFiniteElements.Operators
result, _ = operators.pDtq(upper_bounds, lower_bounds, A, B, C,omega) # Correct call
println(result)

size(result)
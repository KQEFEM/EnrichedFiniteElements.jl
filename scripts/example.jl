using Revise  # If you're using Revise (recommended)
using EnrichedFiniteElements

EnrichedFiniteElements.BasisFunctions.phi(0.5, 0.5) # Test a function
println("Test successful!") # To confirm that the code ran

# Example of shorter alias:
const basis = EnrichedFiniteElements.BasisFunctions
basis.p_matrix(0.5, 0.5)
println("Alias Test successful!") # To confirm that the code ran

grad_matrix_input = [1.0 2.0 3.0; 4.0 5.0 6.0; 7.0 8.0 9.0] # Example matrix
grad_matrix = basis.grads_matrix(grad_matrix_input)
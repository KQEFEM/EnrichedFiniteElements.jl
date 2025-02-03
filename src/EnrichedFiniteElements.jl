module EnrichedFiniteElements

include("BasisFunctions.jl")  # Include BasisFunctions FIRST
using .BasisFunctions       # Use BasisFunctions

include("IntegrationWrappers.jl") # Include IntegrationWrappers SECOND
using .IntegrationWrappers  # Use IntegrationWrappers

include("Operators.jl")      # Include Operators LAST
using .Operators           # Use Operators

end # module EnrichedFiniteElements
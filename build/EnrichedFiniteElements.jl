module EnrichedFiniteElements
export generate_mesh, plot_mesh  # Explicitly export functions

include("BasisFunctions.jl")  # Include BasisFunctions FIRST
using .BasisFunctions       # Use BasisFunctions

include("IntegrationWrappers.jl") # Include IntegrationWrappers SECOND
using .IntegrationWrappers  # Use IntegrationWrappers

include("Operators.jl")      # Include Operators LAST
using .Operators           # Use Operators

include("MatrixCreation.jl")
using .MatrixCreation

include("MeshCreation.jl")
using .MeshCreation

include("EnrichmentCreator.jl")
using .EnrichmentCreator

end  # module EnrichedFiniteElements

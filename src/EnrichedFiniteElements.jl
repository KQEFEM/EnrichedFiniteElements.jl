module EnrichedFiniteElements

export generate_mesh, plot_mesh  # Explicitly export functions

using JuliaFormatter

include("BasisFunctions.jl")  # Include BasisFunctions FIRST
using .BasisFunctions       # Use BasisFunctions

include("IntegrationWrappers.jl") # Include IntegrationWrappers SECOND
using .IntegrationWrappers  # Use IntegrationWrappers

include("Operators.jl")      # Include Operators LAST
using .Operators           # Use Operators

include("TransformationFunctions.jl")
using .TransformationFunctions

include("MatrixCreation.jl")
using .MatrixCreation

include("MeshCreation.jl")
using .MeshCreation

include("EnrichmentCreator.jl")
using .EnrichmentCreator

include("Formulation.jl")
using .Formulation



end  # module EnrichedFiniteElements

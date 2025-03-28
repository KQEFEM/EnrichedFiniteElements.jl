using Revise
using EnrichedFiniteElements
using Test

@testset "EnrichedFiniteElements.jl" begin
    include("BasisFunctionsTests.jl")   # Include the tests for BasisFunctions
    include("EnrichmentCreatorTests.jl") # Include the tests for EnrichmentCreator
    include("OperatorsTests.jl")      # Include the tests for Operatorsend
    include("MatrixCreationTests.jl")
    include("TransformationFunctionsTests.jl")
    include("ForumlationTests.jl")
end

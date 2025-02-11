using EnrichedFiniteElements
using Documenter
using DocMeta

# Set up the doctest environment
DocMeta.setdocmeta!(
    EnrichedFiniteElements,
    :DocTestSetup,
    :(using EnrichedFiniteElements);
    recursive = true,
)

# Build the documentation (Markdown format)
makedocs(
    format = Documenter.Markdown(),
    modules = [EnrichedFiniteElements],
    authors = "Kieran Quaine",
    pages = ["Home" => "index.md"],
)

# Deploy the docs (ensure you have a GITHUB_TOKEN set for this to work)
deploydocs(; repo = "github.com/KQEFEM/EnrichedFiniteElements.jl", devbranch = "main")
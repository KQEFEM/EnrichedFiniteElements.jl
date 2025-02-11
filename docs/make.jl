using EnrichedFiniteElements
using Documenter

# Build the documentation (Markdown format)
makedocs(
    modules = [EnrichedFiniteElements],  # Include your module
    authors = "Kieran Quaine",           # Specify the author
    format = Documenter.Markdown(),      # Use Markdown format
    pages = ["Home" => "index.md"],      # Define the pages
)

# Deploy the docs (ensure you have a GITHUB_TOKEN set for this to work)
deploydocs(
    repo = "github.com/KQEFEM/EnrichedFiniteElements.jl",  # Repository URL
    devbranch = "main",                                    # Default branch
)
using EnrichedFiniteElements
using Documenter

DocMeta.setdocmeta!(
    EnrichedFiniteElements,
    :DocTestSetup,
    :(using EnrichedFiniteElements);
    recursive = true,
)

makedocs(;
    modules = [EnrichedFiniteElements],
    authors = "Kieran Quaine",
    sitename = "EnrichedFiniteElements.jl",
    format = Documenter.HTML(;
        canonical = "https://KQEFEM.github.io/EnrichedFiniteElements.jl",
        edit_link = "main",
        assets = String[],
    ),
    pages = ["Home" => "index.md"],
)

deploydocs(; repo = "github.com/KQEFEM/EnrichedFiniteElements.jl", devbranch = "main")

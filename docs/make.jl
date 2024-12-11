using Documenter
using SplineGrids

makedocs(modules = [SplineGrids],
    sitename = "SplineGrids.jl",
    clean = true,
    linkcheck = true,
    pages = ["index.md",
        "Examples" => [
            "Dimensionality" => "examples_dimensions.md",
            "Derivatives" => "examples_derivatives.md"
        ],
        "Manual" => "manual.md"])

deploydocs(; repo = "github.com/SouthEndMusic/SplineGrids.jl.git")
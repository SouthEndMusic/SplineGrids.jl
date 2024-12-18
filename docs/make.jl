using Documenter
using SplineGrids

makedocs(modules = [SplineGrids],
    sitename = "SplineGrids.jl",
    clean = true,
    linkcheck = true,
    pages = ["index.md",
        "Examples" => [
            "Dimensionality" => "examples_dimensions.md",
            "Derivatives" => "examples_derivatives.md",
            "Fitting" => "examples_linear_fitting.md",
            "Control point derivatives with Enzyme" => "examples_enzyme.md"
        ],
        "Advanced examples" => [
            "Solving a PDE" => "examples_pde.md"
        ],
        "Manual" => "manual.md"])

deploydocs(; repo = "github.com/SouthEndMusic/SplineGrids.jl.git")
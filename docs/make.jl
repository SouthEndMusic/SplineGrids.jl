using Documenter
using SplineGrids

makedocs(modules = [SplineGrids],
    sitename = "SplineGrids.jl",
    clean = true,
    linkcheck = true,
    pagesonly = true,
    linkcheck_ignore = ["https://www.ag.jku.at/pubs/2016gjkmss.pdf"],
    pages = [
        "index.md",
        "Theory" => [
            "Refinement" => "theory_refinement.md",
            "Local refinement (THB-splines)" => "theory_local_refinement.md"
        ],
        "Examples" => [
            "Dimensionality" => "examples_dimensions.md",
            "NURBS" => "examples_nurbs.md",
            "Derivatives" => "examples_derivatives.md",
            "Linear fitting" => "examples_linear_fitting.md",
            "Control point derivatives with Enzyme" => "examples_enzyme.md"
        ],
        "Advanced examples" => [
            "Solving a PDE" => "examples_pde.md",
            "Optimizing a lens surface" => "examples_optics.md"
        ],
        "Manual" => "manual.md"
    ]
)

deploydocs(; repo = "github.com/SouthEndMusic/SplineGrids.jl.git")

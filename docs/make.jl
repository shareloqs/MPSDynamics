using Documenter
using MPSDynamics

makedocs(
    clean   = true,
    doctest = true,
    modules = [MPSDynamics],
    highlightsig = true,
    sitename = "MPSDynamics.jl",
    authors = "Angus Dunnett",
    expandfirst = []
)

deploydocs(
    repo = "github.com/angusdunnett/MPSDynamics.git",
    devurl = "docs"
)

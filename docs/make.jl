using Documenter
using MPSDynamics

makedocs(
    clean   = true,
    doctest = true,
    modules = [MPSDynamics],
    highlightsig = true,
    sitename = "MPSDynamics.jl",
    authors = "Angus Dunnett, Thibaut Lacroix, Brieuc Le Dé",
    expandfirst = []
)

deploydocs(
    repo = "github.com/shareloqs/MPSDynamics.git",
    devurl = "docs"
)

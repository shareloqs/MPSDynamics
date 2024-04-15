using Documenter
using MPSDynamics

makedocs(
    clean   = true,
    doctest = true,
    modules = [MPSDynamics],
    highlightsig = true,
    sitename = "MPSDynamics.jl",
    authors = "Angus Dunnett, Thibaut Lacroix, Brieuc Le Dé, Angela Riva",
    pages = [
        "index.md",
        "user-guide.md",
        "Examples" => ["./examples/sbm.md", "./examples/puredephasing.md"],
        "theory.md",
        "Methods" => "methods.md",
        "dev.md"
    ],
    expandfirst = [],
    remotes = nothing
)

deploydocs(
    repo = "github.com/shareloqs/MPSDynamics.git",
    devurl = "docs",
    devbranch = "doc-writing"
)

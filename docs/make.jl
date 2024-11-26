using Documenter
using MPSDynamics

makedocs(
    clean   = true,
    doctest = true,
    modules = [MPSDynamics],
    highlightsig = true,
    sitename = "MPSDynamics.jl",
    authors = "Angus J. Dunnett, Thibaut Lacroix, Brieuc Le Dé, Angela Riva",
    pages = [
        "index.md",
        "user-guide.md",
        "nutshell.md",
        "Examples" => ["./examples/sbm.md", "./examples/puredephasing.md", "./examples/timedep.md", "./examples/anderson-model.md", "./examples/bath-observables.md", "./examples/protontransfer.md"],
        "convergence.md",
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
    devbranch = "master"
)

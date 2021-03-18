using Documenter
using MPSDynamics

#push!(LOAD_PATH,"../include/")

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
    repo = "git@github.com:angusdunnett/MPSDynamics.git",
    devurl = "docs"
)

# bitbucket = "/home/angus/Documents/Julia/angus-dunnett.bitbucket.io/"
# build = "/home/angus/Documents/Julia/mps2.0/docs/build/"

# Base.Filesystem.cp(joinpath(build, "assets"), joinpath(bitbucket,"assets"), force=true)
# Base.Filesystem.cp(joinpath(build,"index.html"), joinpath(bitbucket,"index.html"), force=true)
# Base.Filesystem.cp(joinpath(build,"search"), joinpath(bitbucket,"search"), force=true)
# Base.Filesystem.cp(joinpath(build,"search_index.js"), joinpath(bitbucket,"search_index.js"), force=true)
# cd(bitbucket)
# pwd()
# run(Cmd(`git commit -a -m "auto-commit"`))
# run(Cmd(`git push`))


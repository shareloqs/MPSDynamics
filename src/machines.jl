struct Machine
    name::String
    exename::String
    nproc::Int
end

function init(mach::Machine)
    pid = addprocs([(mach.name, 1)], exename=mach.exename)
    remotecall_fetch(() -> Pkg.add(PackageSpec(url="https://angus-dunnett@bitbucket.org/angus-dunnett/mpsdynamics.git", rev="master")), only(pid))
    rmprocs(pid)
end


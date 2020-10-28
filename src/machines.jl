struct Machine
    name::String
    exename::String
    nproc::Int
    wdir::String
end
struct LocalMachine
    name::String
end
LocalMachine() = LocalMachine("local")

rmworkers() = rmprocs(workers())

function launch_workers(mach::Machine, nworkers::Int=1)
    pids = addprocs([(mach.name, nworkers)], exename=mach.exename, dir=mach.wdir)
    return pids
end
function launch_workers(machs::Vector{Machine}, nworkers::Int=1)
    pids = Int[]
    for mach in machs
        push!(pids, launch_workers(mach, nworkers)...)
    end
    return pids
end
function launch_workers(f::Function, args...)
    pids = launch_workers(args...)
    try
        ret = f(pids)
    finally
        rmprocs(pids)
    end
    return ret
end

function init_machines(machs::Vector{Machine})
    launch_workers(machs) do pid
        @everywhere [pid...] eval(using Pkg)
        @everywhere [pid...] Pkg.add(PackageSpec(url="https://angus-dunnett@bitbucket.org/angus-dunnett/mpsdynamics.git", rev="master"))
    end
end

function update_machines(machs::Vector{Machine})
    launch_workers(machs) do pid
        @everywhere [pid...] eval(using Pkg)
        @everywhere [pid...] Pkg.update("MPSDynamics")
    end
end


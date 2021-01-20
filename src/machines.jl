abstract type Machine end

struct RemoteMachine <: Machine
    name::String
    exename::String
    wdir::String
end
struct LocalMachine <: Machine
    name::String
end
LocalMachine() = LocalMachine("local")

rmworkers() = rmprocs(workers())

function launch_workers(mach::RemoteMachine, nworkers::Int=1)
    pids = addprocs([(mach.name, nworkers)], exename=mach.exename, dir=mach.wdir, tunnel=true)
    return pids
end
launch_workers(::LocalMachine, nworkers::Int=1) = addprocs(nworkers)
launch_workers(nworkers::Int=1) = launch_workers(::LocalMachine, nworkers)

function launch_workers(machs::Vector{T}, nworkers::Int=1) where T <: Machine
    pids = Int[]
    for mach in machs
        push!(pids, launch_workers(mach, nworkers)...)
    end
    return pids
end
function launch_workers(f::Function, args...)
    pids = launch_workers(args...)
    ret = try
        f(pids)
    finally
        rmprocs(filter(x->x!=1, pids))
    end
    return ret
end

function init_machines(machs::Vector{T}) where T <: Machine
    launch_workers(machs) do pid
        @everywhere pid eval(using Pkg)
        @everywhere pid Pkg.add(PackageSpec(url="https://github.com/angusdunnett/MPSDynamics.git", rev="master"))
    end
end

function update_machines(machs::Vector{T}) where T <: Machine
    launch_workers(machs) do pid
        @everywhere pid eval(using Pkg)
        @everywhere pid Pkg.update("MPSDynamics")
    end
end

anguspc = RemoteMachine("anguspc", "julia", "/home/angus/")
alexpc = RemoteMachine("alexpc", "/home/angus/bin/julia", "/home/angus/")
hp = RemoteMachine("hp", "/home/angus/julia-1.4.2/bin/julia", "/home/angus/")
asusnew = RemoteMachine("asusnew", "/home/angus/julia-1.4.2/bin/julia", "/home/angus/")
asusold = RemoteMachine("asusold", "/home/angus/julia-1.4.2/bin/julia", "/home/angus/")

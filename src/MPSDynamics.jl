module MPSDynamics

using JLD, HDF5, Random, Dates, Plots, Printf, Distributed, LinearAlgebra, DelimitedFiles, KrylovKit, TensorOperations, GraphRecipes, SpecialFunctions, Logging

struct DelayedFunction
    f::Function
    params::Vector
end
evaluate(df::DelayedFunction) = df.f(f.params...)
evaluate(any::Any) = any

struct TensorSim
    dt
    T
    A
    H
    savedir
    params
    obs
    convobs
    savemps
    verbose
    save
    saveplot
    timed
    log
    Dmax
    lightcone
    lightconerad
    lightconethresh
    unid
    name
end

include("config.jl")
include("fundamentals.jl")
include("tensorOps.jl")
include("measure.jl")
include("observables.jl")
include("logiter.jl")
include("machines.jl")
include("treeBasics.jl")
include("treeIterators.jl")
include("treeMeasure.jl")
include("treeTDVP.jl")
include("treeDTDVP.jl")
include("mpsBasics.jl")
include("chainTDVP.jl")
include("chainDMRG.jl")
include("models.jl")

include("runtdvp_dynamic.jl")
include("runtdvp_fixed.jl")
include("runtdvp.jl")
include("convtdvp.jl")

function TensorSim(dt, T, A, H;
                   savedir::String = DEFSAVEDIR,
                   params = [],
                   obs = Observable[],
                   convobs = Observable[],
                   savemps = 0,
                   verbose = false,
                   save = false,
                   saveplot = save,
                   timed = false,
                   log = save,
                   Dmax = throw(error("Dmax must be specified")),
                   lightcone=false,
                   lightconerad=2,
                   lightconethresh=DEFLCTHRESH,
                   unid = randstring(5),
                   name = nothing
                   )
    
    if length(obs)==0
        obs = Observable[]
    end
    if length(convobs)==0
        convobs = Observable[]
    end
    TensorSim(dt,T,A,H,savedir,params,obs,convobs,savemps,verbose,save,saveplot,timed,log,Dmax,lightcone,lightconerad,lightconethresh,unid,name)
end

function runsim(sim::TensorSim, mach::Machine)
    remote = typeof(mach) == RemoteMachine
    remote && update_machines([mach])
    if sim.save || sim.saveplot
        if sim.savedir[end] != '/'
            sim.savedir = string(sim.savedir,"/")
        end
        sim.log || error("the run must be logged if output data is saved") 
        isdir(sim.savedir) || error("save directory $sim.savedir doesn't exist")
    end
    if typeof(sim.Dmax) <: Vector
        convcheck = true
        numDmax = length(sim.Dmax)
    else
        convcheck = false
    end

    if sim.log
        open_log(sim, convcheck, mach)
    end
    errorfile = "$(sim.unid).e"
    
    tstart = now()
    A = dat = nothing
    try
        A, dat = launch_workers(mach) do pid
            print("loading MPSDynamics............")
            @everywhere pid eval(using MPSDynamics)
            println("done")
            A, dat = fetch(@spawnat only(pid) MPSDynamics.runtdvp_fixed!(sim.dt, sim.T, sim.A, sim.H,
                                                                         params=sim.params,
                                                                         obs=sim.obs,
                                                                         convobs=sim.convobs,
                                                                         savemps=sim.savemps,
                                                                         verbose=sim.verbose,
                                                                         timed=sim.timed,
                                                                         Dmax=sim.Dmax,
                                                                         lightcone=sim.lightcone,
                                                                         lightconerad=sim.lightconerad,
                                                                         lightconethresh=sim.lightconethresh,
                                                                         unid=sim.unid
                                                                         ))
            sim.save && save_data(sim.savedir, sim.unid, convcheck, dat["data"], convcheck ? dat["convdata"] : nothing, dat["parameters"])
            convcheck && sim.saveplot && save_plot(sim.savedir, sim.unid, dat["data"]["times"], dat["convdata"], sim.Dmax, sim.convobs)
            return A, dat
        end
    catch e
        sim.log && error_log(sim.savedir, sim.unid)
        showerror(stdout, e, catch_backtrace())                
        println()
        sim.log && open(string(sim.savedir, sim.unid, "/", errorfile), "w+") do io
            showerror(io, e, catch_backtrace())
        end
    finally
        output = length(filter(x-> x!=errorfile && x!="info.txt", readdir(string(sim.savedir, sim.unid)))) > 0
        telapsed = canonicalize(Dates.CompoundPeriod(now() - tstart))
        sim.log && close_log(sim.savedir, sim.unid, output, telapsed)
    end
    return A, dat
end
runsim(sim::TensorSim) = runsim(sim, LocalMachine())

function runsim(sims::Vector, machs::Vector)
    nsims = length(sims)
    f=[]
    launch_workers(nsims) do pids
        @everywhere pids eval(using MPSDynamics)
        for (i, pid) in enumerate(pids)
            push!(f, remotecall(rumsim, pid, sims[i], machs[i]))
        end
    end
    wait.(f)
end

export sz, sx, sy, numb, crea, anih, unitcol, unitrow, unitmat

export chaincoeffs_ohmic, spinbosonmpo, methylbluempo, methylbluempo_correlated, methylbluempo_correlated_nocoupling, methylbluempo_nocoupling

export productstatemps, physdims, randmps, bonddims

export measure, OneSiteObservable, TwoSiteObservable, FockError, errorbar

export TensorSim, runsim, DelayedFunction, evaluate

export Machine, LocalMachine, init_machines, update_machines, launch_workers, rmworkers

export alexpc, hp, asusold, asusnew, anguspc

export VarT, VarX, plot, scatter, loaddat

export printlog, noprintlog

export randtree

export readchaincoeffs

end


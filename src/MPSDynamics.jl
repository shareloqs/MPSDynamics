module MPSDynamics

using JLD, Random, Dates, Plots, Printf, Distributed, LinearAlgebra, DelimitedFiles, KrylovKit, ITensors, TensorOperations, GraphRecipes, SpecialFunctions

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
end

include("config.jl")
include("fundamentals.jl")
include("tensorOps.jl")
include("measure.jl")
include("observables.jl")
include("datiters.jl")
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
                   )
    TensorSim(dt,T,A,H,savedir,params,obs,convobs,savemps,verbose,save,saveplot,timed,log,Dmax,lightcone,lightconerad,lightconethresh,unid)
end

function runsim(sim::TensorSim, mach::Machine)
    update_machines([mach])
    if sim.save || sim.saveplot
        if sim.savedir[end] != '/'
            sim.savedir = string(sim.savedir,"/")
        end
        isdir(sim.savedir) || throw("save directory $sim.savedir doesn't exist")
    end
    if typeof(sim.Dmax) <: Vector
        convcheck = true
        numDmax = length(sim.Dmax)
    else
        convcheck = false
    end
    sim.log && (endpos = open_log(sim, convcheck, mach))
    A, dat = launch_workers(mach) do pid
        tstart = now()
        @everywhere [pid] eval(import MPSDynamics)
        A, dat = fetch(@spawnat pid MPSDynamics.runtdvp_fixed!(sim.dt, sim.T, sim.A, sim.H,
                                                   params=sim.params,
                                                   obs=sim.obs,
                                                   convobs=sim.convobs,
                                                   savemps=sim.savemps,
                                                   verbose=sim.verbose,
                                                   save=false,
                                                   saveplot=false,
                                                   log=false,
                                                   timed=sim.timed,
                                                   Dmax=sim.Dmax,
                                                   lightcone=sim.lightcone,
                                                   lightconerad=sim.lightconerad,
                                                   lightconethresh=sim.lightconethresh,
                                                   unid=sim.unid
                                                   ))
        telapsed = canonicalize(Dates.CompoundPeriod(now() - tstart))
        sim.save && save_data(sim.savedir, dat["data"], dat["convdata"])
        sim.saveplot && save_plot(sim.savedir, dat["convdata"], sim.Dmax, sim.convobs)
        sim.log && close_log(sim.savedir, endpos, telapsed)
        return A, dat
    end
    return A, dat
end
runsim(sim::TensorSim) = runtdvp_fixed!(sim.dt, sim.T, sim.A, sim.H,
                                        params=sim.params,
                                        obs=sim.obs,
                                        savedir=sim.savedir,
                                        convobs=sim.convobs,
                                        savemps=sim.savemps,
                                        verbose=sim.verbose,
                                        save=sim.save,
                                        saveplot=sim.saveplot,
                                        log=sim.log,
                                        timed=sim.timed,
                                        Dmax=sim.Dmax,
                                        lightcone=sim.lightcone,
                                        lightconerad=sim.lightconerad,
                                        lightconethresh=sim.lightconethresh,
                                        unid=sim.unid
                                        )
runsim(sim::TensorSim, mach::LocalMachine) = runsim(sim)
    

export sz, sx, sy, numb, crea, anih, unitcol

export chaincoeffs_ohmic, spinbosonmpo, methylbluempo

export productstatemps, physdims

export measure, OneSiteObservable, TwoSiteObservable

export TensorSim, runsim

export Machine, init_machines, update_machines, launch_workers
    
end

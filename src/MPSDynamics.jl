module MPSDynamics

using JLD, Random, Dates, Plots, Printf, Distributed, LinearAlgebra, DelimitedFiles, KrylovKit, ITensors, TensorOperations, GraphRecipes, SpecialFunctions

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

function run(sims::Vector{TensorSim}, machines::Vector{Machine})
    for mach in machines
        addprocs([(mach.name, mach.nproc)], exename=mach.exename)
    end
    @everywhere eval(using MPSDynamics) # replace with git repo
#    runtdvp_fixed!()
end

f(x)=x-3

export
    sz,
    sx,
    sy,
    numb,
    crea,
    anih,
    unitcol,
    chaincoeffs_ohmic,
    spinbosonmpo,
    methylenebluempo,
    productstatemps,
    TensorSim,
    run,
    runtdvp_fixed!,
    runtdvp_dynamic!,
    physdims,
    measure1siteoperator,
    measure2siteoperator,
    measurempo,
    OneSiteObervable,
    TwoSiteObservable,
    Ntot,
    Nup,
    Ndown,
    SZ,
    SX,
    SY,
    Machine,
    init
end

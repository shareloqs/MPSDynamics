module MPSDynamics

using JLD, HDF5, Random, Dates, Plots, Printf, Distributed, LinearAlgebra, DelimitedFiles, KrylovKit, TensorOperations, GraphRecipes, SpecialFunctions

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
include("chain2TDVP.jl")
include("chainDMRG.jl")
include("models.jl")
include("logging.jl")
include("run_all.jl")
include("run_1TDVP.jl")
include("run_2TDVP.jl")
 
function runsim(dt, tmax, A, H;
                method=:TDVP1,
                machine=LocalMachine(),
                params=[],
                obs=[],
                convobs=[],
                convparams=error("Must specify convergence parameters"),
                save=true,
                saveplot=save,
                savedir="~/MPSDynamics/",
                unid=randstring(5),
                name=nothing,
                kwargs...
    )
    remote = typeof(machine) == RemoteMachine
    remote && update_machines([machine])
    if save || saveplot
        if savedir[end] != '/'
            savedir = string(savedir,"/")
        end
        isdir(savedir) || mkdir(savedir)
        open_log(dt, T, convparams, method, machine, savedir, unid, name, params, obs, convobs, convcheck)
    end
    if typeof(convparams) <: Vector
        convcheck = true
        numconv = length(convparams)
    else
        convcheck = false
    end

    paramdict = Dict([(par[1], par[2]) for par in params]...,
                     [
                         ("dt",dt),
                         ("tmax",tmax),
                         ("method",method),
                         ("convparams",convparams),
                         ("unid",unid),
                         ("name",name)
                     ]
                     )

    errorfile = "$(unid).e"
    
    tstart = now()
    A = dat = nothing
    try
        A, dat = launch_workers(mach) do pid
            print("loading MPSDynamics............")
            @everywhere pid eval(using MPSDynamics)
            println("done")
            A, dat = fetch(@spawnat only(pid) run_all(dt, tmax, A, H;
                                                      method=method,
                                                      obs=obs,
                                                      convobs=convobs,
                                                      convparams=convparams,
                                                      kwargs...))
            save && save_data(savedir, unid, convcheck, dat["data"], convcheck ? dat["convdata"] : nothing, paramdict)
            convcheck && saveplot && save_plot(savedir, unid, dat["data"]["times"], dat["convdata"], convparams, convobs)
            return A, dat
        end
    catch e
        save && error_log(savedir, unid)
        showerror(stdout, e, catch_backtrace())                
        println()
        save && open(string(savedir, unid, "/", errorfile), "w+") do io
            showerror(io, e, catch_backtrace())
        end
    finally
        output = length(filter(x-> x!=errorfile && x!="info.txt", readdir(string(savedir, unid)))) > 0
        telapsed = canonicalize(Dates.CompoundPeriod(now() - tstart))
        save && close_log(savedir, unid, output, telapsed)
    end
    return A, dat
end

export sz, sx, sy, numb, crea, anih, unitcol, unitrow, unitmat

export chaincoeffs_ohmic, spinbosonmpo, methylbluempo, methylbluempo_correlated, methylbluempo_correlated_nocoupling, methylbluempo_nocoupling

export productstatemps, physdims, randmps, bonddims

export measure, OneSiteObservable, TwoSiteObservable, FockError, errorbar

export runsim, run_all

export Machine, RemoteMachine, LocalMachine, init_machines, update_machines, launch_workers, rmworkers

export plot, scatter

export randtree

export readchaincoeffs, h5read

export println, print, show

end


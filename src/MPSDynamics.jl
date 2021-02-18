module MPSDynamics

using JLD, HDF5, Random, Dates, Plots, Printf, Distributed, LinearAlgebra, DelimitedFiles, KrylovKit, TensorOperations, GraphRecipes, SpecialFunctions

include("fundamentals.jl")
include("reshape.jl")
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
include("run_DTDVP.jl")
 
function runsim(dt, tmax, A, H;
                method=:TDVP1,
                machine=LocalMachine(),
                params=[],
                obs=[],
                convobs=[],
                convparams=error("Must specify convergence parameters"),
                save=true,
                plot=save,
                savedir=string(homedir(),"/MPSDynamics/"),
                unid=randstring(5),
                name=nothing,
                kwargs...
                )
    remote = typeof(machine) == RemoteMachine
    remote && update_machines([machine])

    if typeof(convparams) <: Vector && length(convparams) > 1
        convcheck = true
        numconv = length(convparams)
    else
        convcheck = false
        convparams = only(convparams)
    end

    if save || saveplot
        if savedir[end] != '/'
            savedir = string(savedir,"/")
        end
        isdir(savedir) || mkdir(savedir)
        open_log(dt, tmax, convparams, method, machine, savedir, unid, name, params, obs, convobs, convcheck)
    end

    paramdict = Dict([[(par[1], par[2]) for par in params]...,
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
    A0, dat = try
        A0, dat = launch_workers(machine) do pid
            print("loading MPSDynamics............")
            @everywhere pid eval(using MPSDynamics)
            println("done")
            A0, dat = fetch(@spawnat only(pid) run_all(dt, tmax, A, H;
                                                       method=method,
                                                       obs=obs,
                                                       convobs=convobs,
                                                       convparams=convparams,
                                                       kwargs...))
            return A0, dat
        end
        save && save_data(savedir, unid, convcheck, dat["data"], convcheck ? dat["convdata"] : nothing, paramdict)
        plot && save_plot(savedir, convcheck, unid, dat["data"]["times"], convcheck ? dat["convdata"] : dat["data"], convparams, convobs)
        return A0, dat
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
    return A0, dat
end

export sz, sx, sy, numb, crea, anih, unitcol, unitrow, unitmat

export chaincoeffs_ohmic, spinbosonmpo, methylbluempo, methylbluempo_correlated, methylbluempo_correlated_nocoupling, methylbluempo_nocoupling, ibmmpo, methylblue_S1_mpo

export productstatemps, physdims, randmps, bonddims

export measure, measurempo, OneSiteObservable, TwoSiteObservable, FockError, errorbar

export runsim, run_all

export Machine, RemoteMachine, LocalMachine, init_machines, update_machines, launch_workers, rmworkers

export plot, scatter

export randtree

export readchaincoeffs, h5read, load

export println, print, show

end


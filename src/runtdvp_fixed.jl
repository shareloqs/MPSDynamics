function runtdvp_fixed!(dt, T, A, H;
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
                        kwargs...
                        )

    tstart = now()

    obs = union(obs, convobs)

    if save || saveplot
        if savedir[end] != '/'
            savedir = string(savedir,"/")
        end
        isdir(savedir) || throw("save directory $savedir doesn't exist")
    end

    if typeof(Dmax) <: Vector
        convcheck = true
        numDmax = length(Dmax)
    else
        convcheck = false
    end

    endpos = log ? open_log(savedir, dt, T, Dmax, unid, params, obs, convobs, convcheck) : nothing
    
    numsteps = length(collect(0:dt:T))-1
    times = [(i-1)*dt for i=1:numsteps]

    datalist = Any[(par[1], par[2]) for par in params]
    push!(datalist, ("times", times))

    if convcheck
        convdata = Vector{Any}(undef, numDmax-1)
        timed && (convttdvp = Vector{Any}(undef, numDmax-1))
        for (i, D) in enumerate(Dmax[1:end-1])
            @printf("Dmax : %i \n", D)
            data = Vector{Any}(undef, numsteps)
            timed && (ttdvp = Vector{Any}(undef, numsteps))

            F=nothing
            A0=deepcopy(A)
            mpsembed!(A0, D)
            if lightcone
                lc = LightCone(A0, lightconerad, lightconethresh)
            else
                lc = nothing
            end
            for tstep=1:numsteps
                @printf("time = %.3f ", times[tstep])
                lightcone && print(", LCE = $(lc.edge)")
                println()
                data[tstep] = measure(A0, convobs)
                if timed
                    val, t, bytes, gctime, memallocs = @timed tdvp1sweep!(dt, A0, H, F, lc; verbose=verbose, kwargs...)
                    A0, F = val
                    ttdvp[tstep] = t
                else
                    A0, F = tdvp1sweep!(dt, A0, H, F, lc, verbose=verbose, kwargs...)
                end
            end
            convdata[i] = data
            timed && (convttdvp[i] = ttdvp)
        end
        
    end

    D = convcheck ? Dmax[end] : Dmax

    data = Vector{Any}(undef, numsteps)
    timed && (ttdvp = Vector{Any}(undef, numsteps))
    
    @printf("Dmax : %i \n", D)
    F=nothing
    mpsembed!(A, D)
    if lightcone
        lc = LightCone(A, lightconerad, lightconethresh)
    else
        lc = nothing
    end

    for tstep=1:numsteps
        @printf("time = %.3f ", times[tstep])
        lightcone && print(", LCE = $(lc.edge)")
        println()
        data[tstep] = measure(A, obs)
        if timed
            val, t, bytes, gctime, memallocs = @timed tdvp1sweep!(dt, A, H, F, lc; verbose=verbose, kwargs...)
            A, F = val
            ttdvp[tstep] = t
        else
            A, F = tdvp1sweep!(dt, A, H, F, lc; verbose=verbose, kwargs...)
        end
    end

    tdata = Vector{Any}(undef, length(obs))
    for ob=1:length(obs)
        m = data[1][ob]
        for t=2:numsteps
            m = cat(m, data[t][ob]; dims=ndims(obs[ob])+1)
        end
        tdata[ob] = m
    end
    
    push!(datalist, [(ob.name, tdata[i]) for (i,ob) in enumerate(obs)]...)
    timed && push!(datalist, ("tdvptime", ttdvp))

    datadict = Dict(datalist)
    dat = [("data", datadict)]
    if convcheck
        lastprec = map(x->datadict[x.name], convobs)

        #cat causes stack overflow for large numbers of time steps! 
        tconvdata = [[[
            cat([convdata[p][t][ob] for t=1:numsteps]..., dims=ndims(convobs[ob])+1)
            for p=1:numDmax-1]..., lastprec[ob]] for ob=1:length(convobs)]

        convdatalist = Any[(ob.name, tconvdata[i]) for (i,ob) in enumerate(convobs)]
        timed && push!(convdatalist, ("tdvptime", [convttdvp..., ttdvp]))
        push!(convdatalist, ("Dmax", Dmax))
        
        convdatadict = Dict(convdatalist)

        saveplot && save_plot(savedir, unid, convdatadict, Dmax, convobs)

        push!(dat, ("convdata", convdatadict))
    end

    save && save_data(savedir, unid, convcheck, datadict, convdatadict)
    telapsed = canonicalize(Dates.CompoundPeriod(now() - tstart))
    log && close_log(savedir, endpos, telapsed)
    return A, Dict(dat)
end

function open_log(savedir, dt, T, Dmax, unid, params, obs, convobs, convcheck, machine=LocalMachine())
    try
        f = open(string(savedir,"info.txt"))
        close(f)
    catch
        touch(string(savedir,"info.txt"))
    end
    
    endpos = open(string(savedir,"info.txt"), append=true) do f
        writeprintln(f); writeprintln(f)
        writeprintln(f, "start time : $(now())")
        writeprintln(f, "unid : $unid")
        writeprintln(f, "running on : $(machine.name)")
        writeprintln(f, "dt = $dt, tmax = $T")
        writeprintln(f, "parameters : ")
        writeprint(f, "\t")
        for par in params
            writeprint(f, string(par[1], " = ", par[2]), ", ")
        end
        writeprintln(f)
        
        writeprintln(f, "observables : ")
        writeprint(f, "\t")
        for ob in obs
            writeprint(f, ob.name, ", ")
        end
        writeprintln(f)

        if convcheck
            writeprintln(f, "convergence observables : ")
            writeprint(f, "\t")
            for ob in convobs
                writeprint(f, ob.name, ", ")
            end
            writeprintln(f)
        end

        writeprintln(f, "Dmax : $Dmax")
        
        println()
        write(f,"\n >>>>>>>>>>>>>>>>>>>>")
        endpos = position(f)
        write(f,"<<<<<<<<<<<<<<<<<<<<                                                        ")
        return endpos
    end
    return endpos
end
open_log(sim::TensorSim, convcheck, mach=LocalMachine()) = open_log(sim.savedir, sim.dt, sim.T, sim.Dmax, sim.unid, sim.params, sim.obs, sim.convobs, convcheck, mach)

function close_log(savedir, endpos, telapsed)
    open(string(savedir,"info.txt"), "a+") do f
        seek(f, endpos)
        write(f,"run completed at $(now())<<<<<<<<<<<<<<<<<<<<\n")
        write(f,string("run time : ", telapsed, "\n"))
        println(string("run time : ", telapsed))
    end
end

function save_data(savedir, unid, convcheck, datadict, convdatadict)
    jldopen(string(savedir,"dat_",unid,".jld"), "w") do file
        write(file, "data", datadict)
        convcheck && write(file, "convdata", convdatadict)
    end
end

function save_plot(savedir, unid, convdatadict, Dmax, convobs)
    default(size = (800,600), reuse = true, title = unid, legend = true)
    for ob in filter(x->ndims(x)==0, convobs)
        if eltype(convdatadict[ob.name][1]) <: Complex
            plt = plot(hcat(convdatadict[ob.name]...); labels=transpose(Dmax), xlabel="Re($(ob.name))", ylabel="Im($(ob.name))");
        else
            plt = plot(times, convdatadict[ob.name]; labels=transpose(Dmax), xlabel="t", ylabel=ob.name);
        end
        savefig(plt, string(savedir,"convplot_",ob.name,"_",unid,".pdf"));
    end
end

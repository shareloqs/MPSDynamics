function runtdvp_fixed!(dt, T, A, H;
                        params = [],
                        obs = Observable[],
                        convobs = Observable[],
                        savemps = 0,
                        verbose = false,
                        timed = false,
                        Dmax = throw(error("Dmax must be specified")),
                        lightcone=false,
                        lightconerad=2,
                        lightconethresh=DEFLCTHRESH,
                        unid = randstring(5),
                        kwargs...
                        )

    obs = union(obs, convobs)

    numsteps = length(collect(0:dt:T))-1
    times = [(i-1)*dt for i=1:numsteps]

    datalist = Any[(par[1], par[2]) for par in params]
    push!(datalist, ("times", times))

    if typeof(Dmax) <: Vector
        convcheck = true
        numDmax = length(Dmax)
    else
        convcheck = false
    end

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
                @printf("%i/%i, t = %.3f ", tstep, numsteps, times[tstep])
                lightcone && print(", LCE = $(lc.edge)")
                println()
                data[tstep] = measure(A0, convobs; t=times[tstep])
                if timed
                    val, t, bytes, gctime, memallocs = @timed tdvp1sweep!(dt, A0, H, F, lc; verbose=verbose, kwargs...)
                    println("\t","ΔT = ", t)
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
        @printf("%i/%i, t = %.3f ", tstep, numsteps, times[tstep])
        lightcone && print(", LCE = $(lc.edge)")
        println()
        data[tstep] = measure(A, obs; t=times[tstep])
        if timed
            val, t, bytes, gctime, memallocs = @timed tdvp1sweep!(dt, A, H, F, lc; verbose=verbose, kwargs...)
            println("\t","ΔT = ", t)
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
    dat = Any[("data", datadict)]
    if convcheck
        lastprec = map(x->datadict[x.name], convobs)

        #cat causes stack overflow for large numbers of time steps!
        if numDmax != 1
            tconvdata =
                [
                    [
                        [cat([convdata[p][t][ob] for t=1:numsteps]..., dims=ndims(convobs[ob])+1) for p=1:numDmax-1]..., lastprec[ob]
                    ] for ob=1:length(convobs)
                ]
        else
            tconvdata = [[lastprec[ob]] for ob=1:length(convobs)]
        end

        convdatalist = Any[(ob.name, tconvdata[i]) for (i,ob) in enumerate(convobs)]
        timed && push!(convdatalist, ("tdvptime", [convttdvp..., ttdvp]))
        push!(convdatalist, ("Dmax", Dmax))
        
        convdatadict = Dict(convdatalist)
        push!(dat, ("convdata", convdatadict))
    end
    return A, Dict(dat)
end

function open_log(savedir, dt, T, Dmax, unid, params, obs, convobs, convcheck, machine=LocalMachine())
    try
        f = open(string(savedir,"log.txt"))
        close(f)
    catch
        touch(string(savedir,"log.txt"))
    end
    
    open(string(savedir,"log.txt"), append=true) do f
        writeprintln(f, "[$(now())] => RUN <$unid> START")
        writeprintln(f, "\t machine : $(machine.name)")
        writeprintln(f, "\t dt = $dt")
        writeprintln(f, "\t tmax = $T")
        writeprint(f, "\t parameters : ")
        for par in params
            writeprint(f, string(par[1], " = ", par[2]), ", ")
        end
        writeprintln(f)        

        writeprint(f, "\t observables : ")
        for ob in obs
            writeprint(f, ob.name, ", ")
        end
        writeprintln(f)

        if convcheck
            writeprint(f, "\t convergence observables : ")
            for ob in convobs
                writeprint(f, ob.name, ", ")
            end
            writeprintln(f)
        end
        writeprintln(f, "\t Dmax : $Dmax")
        writeprintln(f)
    end
end
open_log(sim::TensorSim, convcheck, mach=LocalMachine()) = open_log(sim.savedir, sim.dt, sim.T, sim.Dmax, sim.unid, sim.params, sim.obs, sim.convobs, convcheck, mach)

function error_log(savedir, unid)
    open(string(savedir,"log.txt"), append=true) do f
        write(f, "[$(now())] => RUN <$unid> ERROR\n")
        write(f, "\t see $(unid)/$(unid).e for details\n")
    end
end

function close_log(savedir, unid, output, telapsed)
    open(string(savedir,"log.txt"), append=true) do f
        writeprintln(f, "[$(now())] => RUN <$unid> END")
        if output
            writeprintln(f, "\t output files produced")
        else
            writeprintln(f, "\t no output files produced")
        end
        writeprintln(f, "\t total run time : $telapsed")
        writeprintln(f)
    end
end

function save_data(savedir, unid, convcheck, datadict, convdatadict)
    jldopen(string(savedir,unid,"/","dat_",unid,".jld"), "w") do file
        write(file, "data", datadict)
        convcheck && write(file, "convdata", convdatadict)
    end
end

function save_plot(savedir, unid, times, convdatadict, Dmax, convobs)
    default(size = (800,600), reuse = true, title = unid, legend = true)
    for ob in filter(x->ndims(x)==0, convobs)
        if eltype(convdatadict[ob.name][1]) <: Complex
            plt = plot(hcat(convdatadict[ob.name]...); labels=transpose(Dmax), xlabel="Re($(ob.name))", ylabel="Im($(ob.name))");
        else
            plt = plot(times, convdatadict[ob.name]; labels=transpose(Dmax), xlabel="t", ylabel=ob.name);
        end
        savefig(plt, string(savedir,unid,"/","convplot_",ob.name,"_",unid,".pdf"));
    end
end

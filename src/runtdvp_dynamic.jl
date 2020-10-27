function runtdvp_dynamic!(dt, T, A, H;
                          savedir::String = DEFSAVEDIR,
                          params = [],
                          obs::Vector{Observable} = Observable[],
                          convobs::Vector{Observable} = Observable[],
                          savemps = 0,
                          verbose = false,
                          save = false,
                          saveplot = save,
                          timed = false,
                          log = save,
                          prec = DEFPREC,
                          Dlim = DEFDLIM,
                          Dplusmax = nothing,
                          error = false,
                          unid = randstring(5),
                          kwargs...
                          )

    tstart = now()

    obs = union(obs, convobs)

    if savedir[end] != '/'
        savedir = string(savedir,"/")
    end

    isdir(savedir) || throw("save directory $savedir doesn't exist")
    
    if typeof(prec) <: Vector
        if length(prec)==1
            prec = prec[1]
            convcheck = false
        else
            convcheck = true
            numprec = length(prec)
        end
    else
        convcheck = false
    end

    if log
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

            writeprintln(f, "prec : $prec")
            
            println()
            write(f,"\n >>>>>>>>>>>>>>>>>>>>")
            endpos = position(f)
            write(f,"<<<<<<<<<<<<<<<<<<<<                                                        ")
            return endpos
        end
    end
    
    numsteps = length(collect(0:dt:T))-1
    times = [(i-1)*dt for i=1:numsteps]

    if convcheck
        convdata = Vector{Any}(undef, numprec-1)
        convdims = Vector{Any}(undef, numprec-1)
        error && (converrs = Vector{Any}(undef, numprec-1))
        timed && (convttdvp = Vector{Any}(undef, numprec-1))
        timed && (convtproj = Vector{Any}(undef, numprec-1))
        for (i, p) in enumerate(prec[1:end-1])
            @printf("prec : %.5e \n", p)
            data = Vector{Any}(undef, numsteps)
            dims = Vector{Any}(undef, numsteps)
            error && (errs = Vector{Any}(undef, numsteps))
            timed && (ttdvp = Vector{Any}(undef, numsteps))
            timed && (tproj = Vector{Any}(undef, numsteps))

            F=nothing
            Afull=nothing
            A0=deepcopy(A)
            for tstep=1:numsteps
                @printf("time = %.3f; Dmax = %i \n", times[tstep], max(bonddimsmps(A0)...))
                A0, Afull, F, info = tdvp1sweep_dynamic!(dt, A0, H, Afull, F;
                                                         obs=convobs,
                                                         prec=p,
                                                         Dlim=Dlim,
                                                         Dplusmax=Dplusmax,
                                                         verbose=verbose,
                                                         timed=timed,
                                                         error=error,
                                                         kwargs...
                                                         )
                data[tstep] = info["obs"]
                dims[tstep] = info["dims"]
                error && (errs[tstep] = info["err"])
                timed && (ttdvp[tstep] = info["t2"] + info["t3"])
                timed && (tproj[tstep] = info["t1"])
            end
            convdata[i] = data
            convdims[i] = hcat(dims...)
            error && (converrs[i] = errs)
            timed && (convttdvp[i] = ttdvp)
            timed && (convtproj[i] = tproj)
        end
    end

    convcheck && (p = prec[end])

    data = Vector{Any}(undef, numsteps)
    dims = Vector{Any}(undef, numsteps)
    efft = Vector{Any}(undef, numsteps)
    error && (errs = Vector{Any}(undef, numsteps))
    timed && (ttdvp = Vector{Any}(undef, numsteps))
    timed && (tproj = Vector{Any}(undef, numsteps))

    @printf("\nprec : %.5e \n\n", p)
    F=nothing
    Afull=nothing
    for tstep=1:numsteps
        @printf("time = %.3f; Dmax = %i \n", times[tstep], max(bonddimsmps(A)...))
        A, Afull, F, info = tdvp1sweep_dynamic!(dt, A, H, Afull, F;
                                                obs=obs,
                                                prec=p,
                                                Dlim=Dlim,
                                                Dplusmax=Dplusmax,
                                                verbose=verbose,
                                                timed=timed,
                                                error=error,
                                                kwargs...
                                                )
        data[tstep] = info["obs"]
        dims[tstep] = info["dims"]
        efft[tstep] = info["effect"]
        error && (errs[tstep] = info["err"])
        timed && (ttdvp[tstep] = info["t2"] + info["t3"])
        timed && (tproj[tstep] = info["t1"])
    end

    tdata = Vector{Any}(undef, length(obs))
    for ob=1:length(obs)
        m = data[1][ob]
        for t=2:numsteps
            m = cat(m, data[t][ob]; dims=ndims(obs[ob])+1)
        end
        tdata[ob] = m
    end

    dims = hcat(dims...)
    datalist = Any[[(ob.name, tdata[i]) for (i,ob) in enumerate(obs)]..., ("bdims", dims), ("effts", efft)]
    
    error && push!(datalist, ("projerr", errs))
    timed && push!(datalist, ("projtime", tproj))
    timed && push!(datalist, ("tdvptime", ttdvp))
    push!(datalist, [(par[1], par[2]) for par in params]...)
    push!(datalist, ("times", times))

    datadict = Dict(datalist)
    if convcheck
        lastprec = map(x->datadict[x.name], convobs)

        #cat causes stack overflow for large numbers of time steps! 
        tconvdata = [[[
            cat([convdata[p][t][ob] for t=1:numsteps]..., dims=ndims(obs[ob])+1)
            for p=1:numprec-1]..., lastprec[ob]] for ob=1:length(convobs)]

        convdatalist = Any[[(ob.name, tconvdata[i]) for (i,ob) in enumerate(convobs)]..., ("bdims", [convdims..., dims])]
        error && push!(convdatalist, ("projerr", [converrs..., errs]))
        timed && push!(convdatalist, ("projtime", [convtproj..., tproj]))
        timed && push!(convdatalist, ("tdvptime", [convttdvp..., ttdvp]))
        push!(convdatalist, ("prec", prec))
        
        convdatadict = Dict(convdatalist)
    end

    if save
        jldopen(string(savedir,"dat_",unid,".jld"), "w") do file
            write(file, "data", datadict)
            convcheck && write(file, "convdata", convdatadict)
        end
    end

    telapsed = canonicalize(Dates.CompoundPeriod(now() - tstart))

    if log
        open(string(savedir,"info.txt"), "a+") do f
            seek(f, endpos)
            write(f,"run completed at $(now())<<<<<<<<<<<<<<<<<<<<\n")
            write(f,string("run time : ", telapsed, "\n"))
            println(string("run time : ", telapsed))
        end
    end

    if convcheck
        return A, datadict, convdatadict
    else
        return A, datadict
    end
end

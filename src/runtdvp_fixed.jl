#Deprecated#
function runtdvp_fixed!(dt, T, A, H;
                        params = [],
                        obs = [],
                        convobs = [],
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

    paramdatalist = Any[(par[1], par[2]) for par in params]
    push!(paramdatalist, ("dt", dt))
    push!(paramdatalist, ("tmax", T))
    datalist = Any[("times", times)]

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

        tconvdata = Vector{Any}(undef, length(convobs))
        for ob=1:length(convobs)
            mp = Vector{Any}(undef, numDmax)
            for p=1:numDmax-1
                m = convdata[p][1][ob]
                for t=2:numsteps
                    m = cat(m, convdata[p][t][ob]; dims=ndims(convobs[ob])+1)
                end
                mp[p] = m
            end
            mp[end] = lastprec[ob]
            tconvdata[ob] = cat(mp...; dims=ndims(convobs[ob])+2)
        end

        convdatalist = Any[(ob.name, tconvdata[i]) for (i,ob) in enumerate(convobs)]
        timed && push!(convdatalist, ("tdvptime", [convttdvp..., ttdvp]))
        push!(convdatalist, ("Dmax", Dmax))
        
        convdatadict = Dict(convdatalist)
        push!(dat, ("convdata", convdatadict))
    end
    push!(dat, ("parameters", Dict(paramdatalist)))
    return A, Dict(dat)
end


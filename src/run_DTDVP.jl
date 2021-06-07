function run_DTDVP(dt, tmax, A, H, prec; obs=[], effects=false, error=false, timed=false, savebonddims=false, Dplusmax=nothing, Dlim=100, kwargs...)
    A0=deepcopy(A)
    data = Dict{String,Any}()

    numsteps = length(collect(0:dt:tmax))-1
    times = [(i-1)*dt for i=1:numsteps+1]
    
    @printf("prec : %.3e \n", prec)

    exp = measure(A0, obs; t=times[1])
    for i=1:length(obs)
        push!(data, obs[i].name => reshape(exp[i], size(exp[i])..., 1))
    end

    bonds = bonddims(A0)
    savebonddims && push!(data, "bonddims" => reshape([bonds...], length(bonds), 1))

    error && (errs = Vector{Float64}(undef, numsteps))
    timed && (ttdvp = Vector{Float64}(undef, numsteps))
    timed && (tproj = Vector{Float64}(undef, numsteps))
    effects && (efft = Vector{Any}(undef, numsteps))

    F=nothing
    Afull=nothing
    for tstep=1:numsteps
        maxbond = max(bonds...)
        @printf("%i/%i, t = %.3f, Dmax = %i \n", tstep, numsteps, times[tstep], maxbond)
        A0, Afull, F, info = tdvp1sweep_dynamic!(dt, A0, H, Afull, F;
                                                 obs=obs,
                                                 prec=prec,
                                                 Dlim=Dlim,
                                                 Dplusmax=Dplusmax,
                                                 timed=timed,
                                                 error=error,
                                                 kwargs...
                                                 )
        
        exp = info["obs"]
        bonds = info["dims"]
        effects && (efft[tstep] = reduce(hcat, info["effect"]))
        error && (errs[tstep] = info["err"])
        timed && (ttdvp[tstep] = info["t2"] + info["t3"])
        timed && (tproj[tstep] = info["t1"])

        if tstep != 1
            for (i, ob) in enumerate(obs)
                data[ob.name] = cat(data[ob.name], exp[i]; dims=ndims(exp[i])+1)
            end
        end
        if savebonddims
            data["bonddims"] = cat(data["bonddims"], bonds, dims=2)
        end
    end
    exp = measure(A0, obs; t=times[end])
    for (i, ob) in enumerate(obs)
        data[ob.name] = cat(data[ob.name], exp[i]; dims=ndims(exp[i])+1)
    end
    
    if effects
        efftarray = efft[1]
        for tstep=2:numsteps
            efftarray = cat(efftarray, efft[tstep], dims=3)
        end
        push!(data, "effects" => efftarray)
    end
    timed && push!(data, "projtime" => tproj)
    timed && push!(data, "tdvptime" => ttdvp)
    error && push!(data, "projerr" => errs)
    push!(data, "times" => times)
    return A0, data
end

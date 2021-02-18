function run_DTDVP(dt, tmax, A, H, prec; obs=[], effects=false, error=false, timed=false, savebonddims=false, Dplusmax=nothing, Dlim=100, kwargs...)
    A0=deepcopy(A)
    data = Dict{String,Any}()

    numsteps = length(collect(0:dt:tmax))-1
    times = [(i-1)*dt for i=1:numsteps+1]
    
    @printf("prec : %.3e \n", prec)

    exp = measure(A0, obs; t=times[1])
    push!(data, [obs[i].name => reshape(exp[i], size(exp[i])..., 1) for i=1:length(obs)]...)

    bonds = bonddims(A0)
    savebonddims && push!(data, "bonddims" => reshape([bonds...], length(bonds), 1))

    error && (errs = Vector{Any}(undef, numsteps))
    timed && (ttdvp = Vector{Any}(undef, numsteps))
    timed && (tproj = Vector{Any}(undef, numsteps))
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
        effects && (efft[tstep] = info["effect"])
        error && (errs[tstep] = info["err"])
        timed && (ttdvp[tstep] = info["t2"] + info["t3"])
        timed && (tproj[tstep] = info["t1"])

        for (i, ob) in enumerate(obs)
            data[ob.name] = cat(data[ob.name], exp[i]; dims=ndims(exp[i])+1)
        end
        if savebonddims
            data["bonddims"] = cat(data["bonddims"], bonds, dims=2)
        end
    end

    effects && push!(data, "effects" => efft)
    timed && push!(data, "projtime" => tproj)
    timed && push!(data, "tdvptime" => ttdvp)
    error && push!(data, "projerr" => errs)
    push!(data, "times" => times)
    return A0, data
end

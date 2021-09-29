function run_A1TDVP(dt, tmax, A, H, prec=10^-4; obs=[], timed = false, savebonddims = true, kwargs...)
    A0=deepcopy(A)
    data = Dict{String,Any}()

    numsteps = length(collect(0:dt:tmax))-1
    times = [(i-1)*dt for i=1:numsteps+1]
    
    @printf("prec : %f \n", prec)

    exp = measure(A0, obs; t=times[1])
    for i=1:length(obs)
        push!(data, obs[i].name => reshape(exp[i], size(exp[i])..., 1))
    end

    timed && (ttdvp = Vector{Float64}(undef, numsteps))
    projerr = Vector{Float64}(undef, numsteps)
    err = Vector{Float64}(undef, numsteps)

    bonds = bonddims(A0)
    savebonddims && push!(data, "bonddims" => reshape([bonds...], length(bonds), 1))

    F=nothing
    Acomp=nothing
    for tstep=1:numsteps
        maxbond = max(bonds...)
        @printf("%i/%i, t = %.3f, Dmax = %i \n", tstep, numsteps, times[tstep], maxbond)
        if timed
            val, t, bytes, gctime, memallocs = @timed tdvpa1sweep!(dt, A0, H, Acomp, F; prec=prec, kwargs...)
            println("\t","ΔT = ", t)
            A0, Acomp, F, info = val
            ttdvp[tstep] = t
        else
            A0, Acomp, F, info = tdvpa1sweep!(dt, A0, H, Acomp, F; prec=prec, kwargs...)
        end
        exp = measure(A0, obs; t=times[tstep])
        for (i, ob) in enumerate(obs)
            data[ob.name] = cat(data[ob.name], exp[i]; dims=ndims(exp[i])+1)
        end
        bonds = info["dims"]
        if savebonddims
            data["bonddims"] = cat(data["bonddims"], bonds, dims=2)
        end
        projerr[tstep] = info["projerr"]
        err[tstep] = info["err"]
    end
    timed && push!(data, "deltat"=>ttdvp)
    push!(data, "times" => times)
    push!(data, "projerr" => projerr)
    push!(data, "err" => err)
    return A0, data
end


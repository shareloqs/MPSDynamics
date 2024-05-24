function run_2TDVP(dt, tmax, A, H, truncerr; obs=[], Dlim=50, savebonddims=false, timed=false, reduceddensity=false, timedep=false, kwargs...)
    A0=deepcopy(A)
    H0=deepcopy(H)
    data = Dict{String,Any}()

    numsteps = length(collect(0:dt:tmax))-1
    times = [(i-1)*dt for i=1:numsteps+1]
    
    @printf("truncerr : %.3e, truncdim : %i \n", truncerr, Dlim)

    exp = measure(A0, obs; t=times[1])
    for i=1:length(obs)
        push!(data, obs[i].name => reshape(exp[i], size(exp[i])..., 1))
    end
    if reduceddensity
        :Nrho in keys(kwargs) ? Nrho = kwargs[:Nrho] : Nrho = 1
        exprho = rhoreduced_1site(A0,Nrho)
        push!(data, "Reduced ρ" => reshape(exprho, size(exprho)..., 1))
    end

    bonds = bonddims(A0)
    savebonddims && push!(data, "bonddims" => reshape([bonds...], length(bonds), 1))

    timed && (ttdvp = Vector{Float64}(undef, numsteps))

    F=nothing
    for tstep=1:numsteps
        maxbond = max(bonds...)
        @printf("%i/%i, t = %.3f, Dmax = %i ", tstep, numsteps, times[tstep], maxbond)
        println()
        if timedep
           Ndrive = kwargs[:Ndrive]
           Htime = kwargs[:Htime]
           H0[Ndrive][1,1,:,:] = H[Ndrive][1,1,:,:] + Htime[tstep][:,:]
        end
        if timed
            val, t, bytes, gctime, memallocs = @timed tdvp2sweep!(dt, A0, H0, F; truncerr=truncerr, truncdim=Dlim, kwargs...)
            println("\t","ΔT = ", t)
            A0, F = val
            ttdvp[tstep] = t
        else
            A0, F = tdvp2sweep!(dt, A0, H0, F; truncerr=truncerr, truncdim=Dlim, kwargs...)
        end
        bonds = bonddims(A0)
        exp = measure(A0, obs; t=times[tstep])
        for (i, ob) in enumerate(obs)
            data[ob.name] = cat(data[ob.name], exp[i], dims=ndims(exp[i])+1)
        end
        if reduceddensity
            exprho = rhoreduced_1site(A0,Nrho)
            data["Reduced ρ"] = cat(data["Reduced ρ"], exprho; dims=ndims(exprho)+1)
        end
        if savebonddims
            data["bonddims"] = cat(data["bonddims"], [bonds...], dims=2)
        end
    end
    timed && push!(data, "deltat"=>ttdvp)
    push!(data, "times" => times)
    return A0, data
end

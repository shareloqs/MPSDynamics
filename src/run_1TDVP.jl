function run_1TDVP(dt, tmax, A, H, Dmax; obs=[], timed=false, reduceddensity=false, timedep=false, kwargs...)
    A0=deepcopy(A)
    H0=deepcopy(H)
    data = Dict{String,Any}()

    numsteps = length(collect(0:dt:tmax))-1
    times = [(i-1)*dt for i=1:numsteps+1]
    
    @printf("Dmax : %i \n", Dmax)

    exp = measure(A0, obs; t=times[1])
    for i=1:length(obs)
        push!(data, obs[i].name => reshape(exp[i], size(exp[i])..., 1))
    end
    if reduceddensity
        :Nrho in keys(kwargs) ? Nrho = kwargs[:Nrho] : Nrho = 1
        exprho = rhoreduced_1site(A0,Nrho)
        push!(data, "Reduced ρ" => reshape(exprho, size(exprho)..., 1))
    end

    timed && (ttdvp = Vector{Float64}(undef, numsteps))

    F=nothing
    mpsembed!(A0, Dmax)
    for tstep=1:numsteps
        @printf("%i/%i, t = %.3f ", tstep, numsteps, times[tstep])
        println()
        if timedep
	   Ndrive = kwargs[:Ndrive]
	   Htime = kwargs[:Htime]
           if length(Ndrive)==1
              size(H0[Ndrive][1,1,:,:])==size(Htime[tstep][:,:]) ? nothing : throw(error("The size of Htime does not match the size of the non-interacting part of H at Ndrive"))
              H0[Ndrive][1,1,:,:] = H[Ndrive][1,1,:,:] + Htime[tstep][:,:]
           else
              for i=1:length(Ndrive)
                 site=Ndrive[i]
                 size(H0[site][end,1,:,:])==size(Htime[site][tstep][:,:]) ? nothing : throw(error("The size of Htime does not match the size of the non-interacting part of H at Ndrive=$site"))
                 H0[site][end,1,:,:] = H[site][end,1,:,:] + Htime[site][tstep][:,:]
              end
           end
        end
        if timed
            val, t, bytes, gctime, memallocs = @timed tdvp1sweep!(dt, A0, H0, F; kwargs...)
            println("\t","ΔT = ", t)
            A0, F = val
            ttdvp[tstep] = t
        else
            A0, F = tdvp1sweep!(dt, A0, H0, F; kwargs...)
        end
        exp = measure(A0, obs; t=times[tstep])
        for (i, ob) in enumerate(obs)
            data[ob.name] = cat(data[ob.name], exp[i]; dims=ndims(exp[i])+1)
        end
        if reduceddensity
            exprho = rhoreduced_1site(A0,Nrho)
            data["Reduced ρ"] = cat(data["Reduced ρ"], exprho; dims=ndims(exprho)+1)
        end
    end
    timed && push!(data, "deltat"=>ttdvp)
    push!(data, "times" => times)
    return A0, data
end

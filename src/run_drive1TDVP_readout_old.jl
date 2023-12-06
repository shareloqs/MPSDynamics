function run_drive1TDVP_readout(dt, tmax, A, H, Dmax; obs=[], Htime, dR, measureparamsChain, timed=false, kwargs...)
    
    A0=deepcopy(A)
    data = Dict{String,Any}()

    H0=deepcopy(H)


    numsteps = length(collect(0:dt:tmax))-1
    times = [(i-1)*dt for i=1:numsteps+1]
    
    @printf("Dmax : %i \n", Dmax)

    exp = measure(A0, obs; t=times[1])
    expChain = measure2siteoperator(A0, measureparamsChain[1], measureparamsChain[2], measureparamsChain[3])


    for i=1:length(obs)
        push!(data, obs[i].name => reshape(exp[i], size(exp[i])..., 1))
    end
    
    push!(data, "Measure Matrix Chain" => reshape(expChain, size(expChain)..., 1))
    timed && (ttdvp = Vector{Float64}(undef, numsteps))

    F=nothing
    mpsembed!(A0, Dmax)
    for tstep=1:numsteps
        
	@printf("%i/%i, t = %.3f ", tstep, numsteps, times[tstep])
        println()
	for i=1:dR
		for j=1:dR
                	H0[2][3,1,i,j] =H[2][3,1,i,j]+Htime[tstep][i,j] #no rwa
			#Drive the cavity for the MPO with chain
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

        expChain = measure2siteoperator(A0, measureparamsChain[1],measureparamsChain[2], measureparamsChain[3])
        data["Measure Matrix Chain"] = cat(data["Measure Matrix Chain"], expChain; dims=ndims(expChain)+1)
	expAtot = A0
        if tstep ==1
                push!(data, "Amplitude Atot" => reshape(expAtot, size(expAtot)..., 1))
        end
	
	data["Amplitude Atot"] = cat(data["Amplitude Atot"], expAtot; dims=ndims(expAtot)+1)


        for (i, ob) in enumerate(obs)
            data[ob.name] = cat(data[ob.name], exp[i]; dims=ndims(exp[i])+1)
        end
    end
    timed && push!(data, "deltat"=>ttdvp)
    push!(data, "times" => times)
    return A0, data
end

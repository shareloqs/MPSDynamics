#Add the include file line in the MPSDynamics.jl
function run_drive1TDVPchain_readout(dt, tmax, A, H, Dmax; obs=[], Htime, dchain, chainparams, measureparamsChain, timed=false, kwargs...)
    
    A0=deepcopy(A)
    data = Dict{String,Any}()

    H0=deepcopy(H)


    numsteps = length(collect(0:dt:tmax))-1
    times = [(i-1)*dt for i=1:numsteps+1]
    
    @printf("Dmax : %i \n", Dmax)

    exp = measure(A0, obs; t=times[1])
    expChain = measure2siteoperator(A0, measureparamsChain[1], measureparamsChain[2], measureparamsChain[3])
    exprho = rhoreduced_proton2chains(A0, 2) #2 to target R (cavity) site, done for the Wigner function


    for i=1:length(obs)
        push!(data, obs[i].name => reshape(exp[i], size(exp[i])..., 1))
    end

    push!(data, "Measure Matrix Chain" => reshape(expChain, size(expChain)..., 1))
    push!(data, "Measure Reduced ρ" => reshape(exprho, size(exprho)..., 1))
    timed && (ttdvp = Vector{Float64}(undef, numsteps))

    F=nothing
    mpsembed!(A0, Dmax)
    tchain = chainparams[2]
    for tstep=1:numsteps
        
	@printf("%i/%i, t = %.3f ", tstep, numsteps, times[tstep])
        println()
	for i=1:dchain
		for j=1:dchain
			for k=1:(length(H)-2)
                		H0[k+2][end,1,i,j] =H[k+2][end,1,i,j]+Htime[k,tstep,i,j] #if no rwa
			#Drive the chain for a particular star frequency
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

        expChain = measure2siteoperator(A0, measureparamsChain[1],measureparamsChain[2], measureparamsChain[3])
        data["Measure Matrix Chain"] = cat(data["Measure Matrix Chain"], expChain; dims=ndims(expChain)+1)
	expAtot = A0
        if tstep ==1
                push!(data, "Amplitude Atot" => reshape(expAtot, size(expAtot)..., 1))
        end
	
	data["Amplitude Atot"] = cat(data["Amplitude Atot"], expAtot; dims=ndims(expAtot)+1)
        exprho = rhoreduced_proton2chains(A0, 2)
        data["Measure Reduced ρ"] = cat(data["Measure Reduced ρ"], exprho; dims=ndims(exprho)+1)


        for (i, ob) in enumerate(obs)
            data[ob.name] = cat(data[ob.name], exp[i]; dims=ndims(exp[i])+1)
        end
    end
    timed && push!(data, "deltat"=>ttdvp)
    push!(data, "times" => times)
    return A0, data
end

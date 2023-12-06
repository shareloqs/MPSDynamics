function run_drive1TDVP(dt, tmax, A, H, Dmax; obs=[], Htime, timed=false, kwargs...)
    
    A0=deepcopy(A)
    data = Dict{String,Any}()

    H0=deepcopy(H)


    numsteps = length(collect(0:dt:tmax))-1
    times = [(i-1)*dt for i=1:numsteps+1]
    
    @printf("Dmax : %i \n", Dmax)

    exp = measure(A0, obs; t=times[1])
    for i=1:length(obs)
        push!(data, obs[i].name => reshape(exp[i], size(exp[i])..., 1))
    end

    timed && (ttdvp = Vector{Float64}(undef, numsteps))

    F=nothing
    mpsembed!(A0, Dmax)
    for tstep=1:numsteps
        for i=1:2
		for j=1:2
        		H0[1][:,:,i,j] = [H[1][:,:,i,j][1]+Htime[tstep][i,j],H[1][:,:,i,j][2],H[1][:,:,i,j][3]] 
		end
	end
        @printf("%i/%i, t = %.3f ", tstep, numsteps, times[tstep])
        println()
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
    end
    timed && push!(data, "deltat"=>ttdvp)
    push!(data, "times" => times)
    return A0, data
end

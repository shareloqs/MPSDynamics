function run_1TDVP(dt, T, A, H, Dmax; obs=[], verbose=false, timed=false)

    numsteps = length(collect(0:dt:T))-1
    times = [(i-1)*dt for i=1:numsteps+1]
    
    @printf("Dmax : %i \n", D)

    exp = measure(A0, obs; t=times[1])
    data = Dict([obs[i].name => exp[i] for i=1:length(obs)])

    timed && (ttdvp = Vector{Float64}(undef, numsteps))

    F=nothing
    A0=deepcopy(A)
    mpsembed!(A0, Dmax)
    for tstep=1:numsteps
        @printf("%i/%i, t = %.3f ", tstep, numsteps, times[tstep])
        println()
        if timed
            val, t, bytes, gctime, memallocs = @timed tdvp1sweep!(dt, A0, H, F; verbose=verbose, kwargs...)
            println("\t","ΔT = ", t)
            A0, F = val
            ttdvp[tstep] = t
        else
            A0, F = tdvp1sweep!(dt, A0, H, F, verbose=verbose, kwargs...)
        end
        exp = measure(A0, obs; t=times[tstep])
        for (i, ob) in enumerate(obs)
            data[ob.name] = cat(data[ob.name], exp[i], dims=ndims(exp[i])+1)
        end
    end
    timed && push!(data, "deltat"=>ttdvp)
    push!(data, "times" => times)
    return A0, data
end

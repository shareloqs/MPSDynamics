
"Old version of run_drive!"
# function run_driveDTDVP_readout(dt, tmax, A, H, prec; obs=[], Htime, dR, drivecav=true, effects=false, error=false, timed=false, savebonddims=false, Dplusmax=nothing, Dlim=50, kwargs...)

"New version with option to only compute the reduced density matrix at some timesteps (input param)"
function run_driveDTDVP_readout(dt, tmax, A, H, prec; 
                                time_frames_rho,    # to specify the time frames at which to compute the reduce ρ
                                obs=[], 
                                Htime, 
                                dR, 
                                drivecav=true, 
                                effects=false, 
                                error=false, 
                                timed=false, 
                                savebonddims=false, 
                                Dplusmax=nothing, 
                                Dlim=50, 
                                kwargs...)
    A0 = deepcopy(A)
    data = Dict{String,Any}()

    H0 = deepcopy(H)

    numsteps = length(collect(0:dt:tmax))-1
    times = [(i-1) * dt for i=1:numsteps + 1]
    
    @printf("prec : %.3e \n", prec)

    exp = measure(A0, obs; t = times[1])
    #expChain = measure2siteoperator(A0, measureparamsChain[1], measureparamsChain[2], measureparamsChain[3])
    exprho = rhoreduced_proton2chains(A0, 2) #2 to target R (cavity) site, done for the Wigner function

    for i=1:length(obs)
        push!(data, obs[i].name => reshape(exp[i], size(exp[i])..., 1))
    end
    #push!(data, "Measure Matrix Chain" => reshape(expChain, size(expChain)..., 1))
    # this is to obtain the reduced density matrix of the cavity
    push!(data, "ρ_red" => reshape(exprho, size(exprho)..., 1))

    bonds = bonddims(A0)
    savebonddims && push!(data, "bonddims" => reshape([bonds...], length(bonds), 1))

    error && (errs = Vector{Float64}(undef, numsteps))
    timed && (ttdvp = Vector{Float64}(undef, numsteps))
    timed && (tproj = Vector{Float64}(undef, numsteps))
    effects && (efft = Vector{Any}(undef, numsteps))

    F = nothing
    Afull = nothing
    for tstep=1:numsteps
        maxbond = max(bonds...)
        
        @printf("%i/%i, t = %.3f, Dmax = %i \n", tstep, numsteps, times[tstep], maxbond)
        if drivecav==true
            for i=1:dR
                for j=1:dR
                        H0[2][3,1,i,j] =H[2][3,1,i,j]+Htime[tstep][i,j] #Drive the cavity 
                end
            end
        else
            for i=1:2
                for j=1:2
                        H0[1][1,1,i,j] =H[1][1,1,i,j]+Htime[tstep][i,j] #Drive the TLS 
                end
            end
        end


        A0, Afull, F, info = tdvp1sweep_dynamic!(dt, A0, H0, Afull, F;
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

        #expChain = measure2siteoperator(A0, measureparamsChain[1],measureparamsChain[2], measureparamsChain[3])
        #data["Measure Matrix Chain"] = cat(data["Measure Matrix Chain"], expChain; dims=ndims(expChain)+1)
        #expAtot = A0
        #if tstep ==1
        #        push!(data, "Amplitude Atot" => reshape(expAtot, size(expAtot)..., 1))
        #end

        #data["Amplitude Atot"] = cat(data["Amplitude Atot"], expAtot; dims=ndims(expAtot)+1)
        
        if tstep*dt ∈ time_frames_rho
            exprho = rhoreduced_proton2chains(A0, 2)
            data["ρ_red"] = cat(data["ρ_red"], exprho; dims=ndims(exprho)+1)
        end
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


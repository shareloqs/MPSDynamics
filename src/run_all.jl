function run_all(dt, tmax, A, H;
                 method=:TDVP1,
                 obs=[],
                 convobs=[],
                 convparams=error("Must specify convergence parameters"),
                 kwargs...
                 )

    obs = union(obs, convobs)

    if typeof(convparams) <: Vector
        convcheck = true
        numconv = length(convparams)
    else
        convcheck = false
    end

    convdat=nothing
    if convcheck
        for (i, cps) in enumerate(convparams[1:end-1])
            if method==:TDVP1
                B, dat = run_1TDVP(dt, tmax, A, H, cps...; obs=convobs, kwargs...)
            elseif method==:TDVP2
                B, dat = run_2TDVP(dt, tmax, A, H, cps...; obs=convobs, kwargs...)
            elseif method==:TDVP1LC
                B, dat = run_1TDVP_LC(dt, tmax, A, H, cps...; obs=convobs, kwargs...)
            elseif method==:DTDVP
                B, dat = run_DTDVP(dt, tmax, A, H, cps...; obs=convobs, kwargs...)
            elseif method==:driveDTDVP_readout
                B, dat = run_driveDTDVP_readout(dt, tmax, A, H, cps...; obs=convobs, kwargs...)
	    elseif method==:A1TDVP
                B, dat = run_A1TDVP(dt, tmax, A, H, cps...; obs=convobs, kwargs...)
	    elseif method==:driveTDVP1
		B, dat = run_drive1TDVP(dt, tmax, A, H, cps...; obs=convobs, kwargs...)
	    elseif method==:driveTDVP1_readout
		B, dat = run_drive1TDVP_readout(dt, tmax, A, H, cps...; obs=convobs, kwargs...)
	    elseif method==:driveTDVP1chain_readout
		B, dat = run_drive1TDVPchain_readout(dt, tmax, A, H, cps...; obs=convobs, kwargs...)
            elseif method==:hTDVP
                error("method $method not recognised")
                #B, dat = run_hTDVP(dt, tmax, A, H, cps...; obs=convobs, kwargs...)
            else
                error("method $method not recognised")
            end
            if i==1
                convdat = Dict([(item.first, reshape(item.second, size(item.second)..., 1)) for item in dat])
            else
                for item in dat
                    convdat[item.first] = cat(convdat[item.first], item.second, dims=ndims(item.second)+1)
                end
            end
        end
    end

    cps = convcheck ? convparams[end] : convparams
    if method==:TDVP1
        B, dat = run_1TDVP(dt, tmax, A, H, cps...; obs=obs, kwargs...)
    elseif method==:TDVP2
        B, dat = run_2TDVP(dt, tmax, A, H, cps...; obs=obs, kwargs...)
    elseif method==:TDVP1LC
        B, dat = run_1TDVP_LC(dt, tmax, A, H, cps...; obs=obs, kwargs...)
    elseif method==:DTDVP
        B, dat = run_DTDVP(dt, tmax, A, H, cps...; obs=obs, kwargs...)
    elseif method==:driveDTDVP_readout
        B, dat = run_driveDTDVP_readout(dt, tmax, A, H, cps...; obs=obs, kwargs...)
    elseif method==:A1TDVP
        B, dat = run_A1TDVP(dt, tmax, A, H, cps...; obs=obs, kwargs...)
    elseif method==:driveTDVP1
        B, dat = run_drive1TDVP(dt, tmax, A, H, cps...; obs=convobs, kwargs...)
    elseif method==:driveTDVP1_readout
        B, dat = run_drive1TDVP_readout(dt, tmax, A, H, cps...; obs=convobs, kwargs...)
    elseif method==:driveTDVP1chain_readout
        B, dat = run_drive1TDVPchain_readout(dt, tmax, A, H, cps...; obs=convobs, kwargs...)

    elseif method==:hTDVP
        error("method $method not recognised")
        #B, dat = run_hTDVP(dt, tmax, A, H, cps...; obs=obs, kwargs...)
    else
        error("method $method not recognised")
    end
    if convcheck
        for item in keys(convdat)
            convdat[item] = cat(convdat[item], dat[item], dims=ndims(dat[item])+1)
        end
    end

    if convcheck
        data = Dict(["data"=>dat, "convdata"=>convdat])
    else
        data = Dict(["data"=>dat])
    end
    return B, data
end

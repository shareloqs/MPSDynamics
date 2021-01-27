function run_all(dt, T, A, H;
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

    if convcheck
        for (i, cps) in enumerate(convparams[1:end-1])
            if method==:TDVP1
                B, dat = run_1TDVP(dt, T, A, H, cps...; obs=convobs, kwargs...)
            elseif method==:TDVP2
                B, dat = run_2TDVP(dt, T, A, H, cps...; obs=convobs, kwargs...)
            elseif method==:TDVP1LC
                B, dat = run_1TDVP_LC(dt, T, A, H, cps...; obs=convobs, kwargs...)
            elseif method==:DTDVP
                B, dat = run_DTDVP(dt, T, A, H, cps...; obs=convobs, kwargs...)
            elseif method==:hTDVP
                B, dat = run_hTDVP(dt, T, A, H, cps...; obs=convobs, kwargs...)
            else
                error("method $method not recognised")
            end
            if i==1
                convdat = deepcopy(dat)
            else
                for item in keys(dat)
                    convdat[item] = cat(convdat[item], dat[item], dims=ndims(dat[item])+1)
                end
            end
        end
    end

    cps = convcheck ? convparams[end] : convparams
    if method==:TDVP1
        B, dat = run_1TDVP(dt, T, A, H, cps...; obs=obs, kwargs...)
    elseif method==:TDVP2
        B, dat = run_2TDVP(dt, T, A, H, cps...; obs=obs, kwargs...)
    elseif method==:TDVP1LC
        B, dat = run_1TDVP_LC(dt, T, A, H, cps...; obs=obs, kwargs...)
    elseif method==:DTDVP
        B, dat = run_DTDVP(dt, T, A, H, cps...; obs=obs, kwargs...)
    elseif method==:hTDVP
        B, dat = run_hTDVP(dt, T, A, H, cps...; obs=obs, kwargs...)
    else
        error("method $method not recognised")
    end
    if convcheck
        for item in keys(convdat)
            convdat[item] = cat(convdat[item], dat[item], dims=ndims(dat[item])+1)
        end
    end

    data = Dict(["data" => dat])
    convcheck && push!(data, "convdata" => convdat)
    return B, data
end

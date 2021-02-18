function runtdvp!(dt, T, B, H, dir::String; params=[], obs=[], saveevery=nothing, verbose=false, save=false, timed=false, unid=nothing, log=true, lightcone=false)

    if !(save)
        return runtdvp!(dt, T, B, H, obs=obs, verbose=verbose, timed=timed, lightcone=lightcone)
    end
    
    if lightcone
        lc = LightCone(B)
    else
        lc = nothing
    end

    if dir[end] != '/'
        dir = string(dir,"/")
    end

    mkpath(dir)

    numobs = length(obs)
    numparams = length(params)
    if unid==nothing
        unid = randstring(5)
    end
    F = nothing
    t = 0
    numsteps = Int64(round(T/dt))

    keylist = [obs[i][1] for i=1:numobs]

    if saveevery == nothing
        saveevery = numsteps*2 # will save only initial and final mps
    end
    
    try
        f = open(string(dir,"info.txt"))
        close(f)
    catch
        touch(string(dir,"info.txt"))
    end

    if log
        endpos = open(string(dir,"info.txt"), append=true) do f
            write(f,"\n\n")
            println("start time : $(now())")
            write(f,"start time : $(now()) \n")
            println("unid : $unid")
            write(f,"unid : $unid \n")
            for i in 1:numparams
                println(string(params[i][1], " : ",params[i][2]))
                write(f, string(params[i][1], " : ",params[i][2], " \n"))
            end
            print("observables : ")
            write(f,"observables : ")
            for i in 1:numobs
                print(keylist[i])
                print(", ")
                write(f, string(keylist[i], ", "))
            end
            print("\n")
            write(f, "\n")
            
            println()
            write(f,"\n >>>>>>>>>>>>>>>>>>>>")
            endpos = position(f)
            write(f,"<<<<<<<<<<<<<<<<<<<<                                                        ")
            return endpos
        end
    end
    
    res = Array{Any}(undef, numobs, numsteps+1)
    times = Vector{Any}()
    savesteps = Vector{Any}()

    tstart=now()
    jldopen(string(dir,"dat_",unid,".jld"), "w") do file
        for i in 1:length(params)
            write(file, params[i][1], params[i][2])
        end
        write(file, "obs", keylist)

        if lightcone
            println("time = 0, LCE = $(lc.edge)")
        else
            println("time = 0")
        end

        for tstep in 1:numsteps
            for i in 1:numobs
                res[i,tstep] = obs[i][2](B)
            end
            if tstep%saveevery == 1
                write(file, "tmps$tstep", B)
                push!(savesteps, tstep)
            end
            push!(times, t)

            B, F = timed ? tdvp1sweep_timed!(dt, B, H, F, lc, verbose=verbose) : tdvp1sweep!(dt, B, H, F, lc, verbose=verbose)
            t += dt
            if lightcone
                println("time = $(round(t;digits=2)), LCE = $(lc.edge)")
            else
                println("time = $(round(t;digits=2))")
            end
        end

        for i in 1:numobs
            res[i,numsteps+1] = obs[i][2](B)
        end
        push!(times, t)

        for i in 1:numobs
            write(file, obs[i][1], res[i,:])
        end
        write(file,"mpsfinal", B)
        write(file,"times",times)
        write(file,"savesteps",savesteps)
    end
    if log
        open(string(dir,"info.txt"), "a+") do f
            seek(f, endpos)
            write(f,"run completed at $(now())<<<<<<<<<<<<<<<<<<<<\n")
            write(f,string("run time : ",now()-tstart,"\n"))
            println(string("run time : ",now()-tstart))
        end
    end
    
    return B, Dict([(obs[i][1],res[i,:]) for i in 1:numobs]), unid
end
runtdvp_timed!(dt, T, B, H, dir::String; params=[], obs=[], saveevery=nothing, verbose=false, save=false, timed=false, unid=nothing, log=true, lightcone=false) = @time runtdvp!(dt, T, B, H, dir; params=params, obs=obs, saveevery=saveevery, verbose=verbose, save=save, timed=timed, unid=unid, log=log, lightcone=lightcone)

function runtdvp!(dt, T, B, H; obs=[], verbose=false, timed=false, lightcone=false)

    if lightcone
        lc = LightCone(B)
    else
        lc = nothing
    end
    
    numobs = length(obs)
    F = nothing
    t = 0
    numsteps = Int64(round(T/dt))
    res = Array{Any}(undef, numobs, numsteps+1)
    
    if lightcone
        println("time = 0, LCE = $(lc.edge)")
    else
        println("time = 0")
    end
    
    for tstep in 1:numsteps
        for i in 1:numobs
            res[i,tstep] = obs[i][2](B)
        end

        B, F = timed ? tdvp1sweep_timed!(dt, B, H, F, lc; verbose=verbose) : tdvp1sweep!(dt, B, H, F, lc; verbose=verbose)
        t += dt
        if lightcone
            println("time = $(round(t;digits=2)), LCE = $(lc.edge)")
        else
            println("time = $(round(t;digits=2))")
        end
    end

    for i in 1:numobs
        res[i,numsteps+1] = obs[i][2](B)
    end
    return B, Dict([(obs[i][1],res[i,:]) for i in 1:numobs])
end
runtdvp_timed!(dt, T, B, H; obs=[], verbose=false, timed=false, lightcone=false) = @time runtdvp!(dt, T, B, H; obs=obs, verbose=verbose, timed=timed, lightcone=lightcone)

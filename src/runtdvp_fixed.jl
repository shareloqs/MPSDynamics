#Deprecated#
function runtdvp_fixed!(dt, T, A, H;
                        params = [],
                        obs = [],
                        convobs = [],
                        savemps = 0,
                        verbose = false,
                        timed = false,
                        Dmax = throw(error("Dmax must be specified")),
                        lightcone=false,
                        lightconerad=2,
                        lightconethresh=DEFLCTHRESH,
                        unid = randstring(5),
                        kwargs...
                        )
    obs = union(obs, convobs)

    numsteps = length(collect(0:dt:T))-1
    times = [(i-1)*dt for i=1:numsteps]

    paramdatalist = Any[(par[1], par[2]) for par in params]
    push!(paramdatalist, ("dt", dt))
    push!(paramdatalist, ("tmax", T))
    datalist = Any[("times", times)]

    if typeof(Dmax) <: Vector
        convcheck = true
        numDmax = length(Dmax)
    else
        convcheck = false
    end

    if convcheck
        convdata = Vector{Any}(undef, numDmax-1)
        timed && (convttdvp = Vector{Any}(undef, numDmax-1))
        for (i, D) in enumerate(Dmax[1:end-1])
            @printf("Dmax : %i \n", D)
            data = Vector{Any}(undef, numsteps)
            timed && (ttdvp = Vector{Any}(undef, numsteps))

            F=nothing
            A0=deepcopy(A)
            mpsembed!(A0, D)
            if lightcone
                lc = LightCone(A0, lightconerad, lightconethresh)
            else
                lc = nothing
            end
            for tstep=1:numsteps
                @printf("%i/%i, t = %.3f ", tstep, numsteps, times[tstep])
                lightcone && print(", LCE = $(lc.edge)")
                println()
                data[tstep] = measure(A0, convobs; t=times[tstep])
                if timed
                    val, t, bytes, gctime, memallocs = @timed tdvp1sweep!(dt, A0, H, F, lc; verbose=verbose, kwargs...)
                    println("\t","ΔT = ", t)
                    A0, F = val
                    ttdvp[tstep] = t
                else
                    A0, F = tdvp1sweep!(dt, A0, H, F, lc, verbose=verbose, kwargs...)
                end
            end
            convdata[i] = data
            timed && (convttdvp[i] = ttdvp)
        end
        
    end

    D = convcheck ? Dmax[end] : Dmax

    data = Vector{Any}(undef, numsteps)
    timed && (ttdvp = Vector{Any}(undef, numsteps))
    
    @printf("Dmax : %i \n", D)
    F=nothing
    mpsembed!(A, D)
    if lightcone
        lc = LightCone(A, lightconerad, lightconethresh)
    else
        lc = nothing
    end

    for tstep=1:numsteps
        @printf("%i/%i, t = %.3f ", tstep, numsteps, times[tstep])
        lightcone && print(", LCE = $(lc.edge)")
        println()
        data[tstep] = measure(A, obs; t=times[tstep])
        if timed
            val, t, bytes, gctime, memallocs = @timed tdvp1sweep!(dt, A, H, F, lc; verbose=verbose, kwargs...)
            println("\t","ΔT = ", t)
            A, F = val
            ttdvp[tstep] = t
        else
            A, F = tdvp1sweep!(dt, A, H, F, lc; verbose=verbose, kwargs...)
        end
    end
    
    tdata = Vector{Any}(undef, length(obs))
    for ob=1:length(obs)
        m = data[1][ob]
        for t=2:numsteps
            m = cat(m, data[t][ob]; dims=ndims(obs[ob])+1)
        end
        tdata[ob] = m
    end
    
    push!(datalist, [(ob.name, tdata[i]) for (i,ob) in enumerate(obs)]...)
    timed && push!(datalist, ("tdvptime", ttdvp))

    datadict = Dict(datalist)
    dat = Any[("data", datadict)]
    if convcheck
        lastprec = map(x->datadict[x.name], convobs)

        tconvdata = Vector{Any}(undef, length(convobs))
        for ob=1:length(convobs)
            mp = Vector{Any}(undef, numDmax)
            for p=1:numDmax-1
                m = convdata[p][1][ob]
                for t=2:numsteps
                    m = cat(m, convdata[p][t][ob]; dims=ndims(convobs[ob])+1)
                end
                mp[p] = m
            end
            mp[end] = lastprec[ob]
            tconvdata[ob] = cat(mp...; dims=ndims(convobs[ob])+2)
        end

        convdatalist = Any[(ob.name, tconvdata[i]) for (i,ob) in enumerate(convobs)]
        timed && push!(convdatalist, ("tdvptime", [convttdvp..., ttdvp]))
        push!(convdatalist, ("Dmax", Dmax))
        
        convdatadict = Dict(convdatalist)
        push!(dat, ("convdata", convdatadict))
    end
    push!(dat, ("parameters", Dict(paramdatalist)))
    return A, Dict(dat)
end

function open_log(savedir, lc, name, dt, T, Dmax, unid, params, obs, convobs, convcheck, machine=LocalMachine())
    mkdir(string(savedir, unid))
    logname = string(savedir,"log-",gethostname(),".txt")
    try
        f = open(logname)
        close(f)
    catch
        touch(logname)
    end

    named = typeof(name) <: String
    open(string(savedir,"$unid/","info.txt"), "w") do f0
        open(logname, append=true) do f
            writeprintln(f, "[$(now())] => RUN <$unid> START")
            named && writeprintln([f,f0], "\t name : $name")
            writeprintln([f,f0], "\t machine : $(machine.name)")
            writeprintln([f,f0], "\t dt = $dt")
            writeprintln([f,f0], "\t tmax = $T")
            writeprintln([f,f0], "\t lightcone : $lc")
            writeprint([f,f0], "\t parameters : ")
            for par in params
                writeprint([f,f0], string(par[1], " = ", par[2]), ", ")
            end
            writeprintln([f,f0])        

            writeprint([f,f0], "\t observables : ")
            for ob in obs
                writeprint([f,f0], ob.name, ", ")
            end
            writeprintln([f,f0])

            if convcheck
                writeprint([f,f0], "\t convergence observables : ")
                for ob in convobs
                    writeprint([f,f0], ob.name, ", ")
                end
                writeprintln([f,f0])
            end
            writeprintln([f,f0], "\t Dmax : $Dmax")
            writeprintln([f,f0])
        end
    end
end
open_log(sim::TensorSim, convcheck, mach=LocalMachine()) = open_log(sim.savedir, sim.lightcone, sim.name, sim.dt, sim.T, sim.Dmax, sim.unid, sim.params, sim.obs, sim.convobs, convcheck, mach)

function error_log(savedir, unid)
    logname = string(savedir,"log-",gethostname(),".txt")
    open(logname, append=true) do f
        write(f, "[$(now())] => RUN <$unid> ERROR\n")
        write(f, "\t see $(unid)/$(unid).e for details\n")
    end
end

function close_log(savedir, unid, output, telapsed)
    logname = string(savedir,"log-",gethostname(),".txt")
    open(logname, append=true) do f
        writeprintln(f, "[$(now())] => RUN <$unid> END")
        if output
            writeprintln(f, "\t output files produced")
        else
            writeprintln(f, "\t no output files produced")
        end
        writeprintln(f, "\t total run time : $telapsed")
        writeprintln(f)
    end
end

function save_data(savedir, unid, convcheck, datadict, convdatadict, paramdatadict)
    jldopen(string(savedir,unid,"/","dat_",unid,".jld"), "w") do file
        g1 = create_group(file, "data")
        for el in datadict
            g1[el.first] = el.second
            if eltype(el.second) <: Complex
                g1[string(el.first,"-re")] = real.(el.second)
                g1[string(el.first,"-im")] = imag.(el.second)
            end
        end
        if convcheck
            g2 = create_group(file, "convdata")
            for el in convdatadict
                g2[el.first] = el.second
                if eltype(el.second) <: Complex
                    g2[string(el.first,"-re")] = real.(el.second)
                    g2[string(el.first,"-im")] = imag.(el.second)
                end
            end
        end
        g3 = create_group(file, "parameters")
        for el in paramdatadict
            g3[el.first] = el.second
            if eltype(el.second) <: Complex
                g3[string(el.first,"-re")] = real.(el.second)
                g3[string(el.first,"-im")] = imag.(el.second)
            end
        end
    end
end

function save_plot(savedir, unid, times, convdatadict, Dmax, convobs)
    default(size = (800,600), reuse = true, title = unid, legend = true)
    for ob in filter(x->ndims(x)==0, convobs)
        if eltype(convdatadict[ob.name]) <: Complex
            plt = plot(convdatadict[ob.name]; labels=transpose(Dmax), xlabel="Re($(ob.name))", ylabel="Im($(ob.name))");
        else
            plt = plot(times, convdatadict[ob.name]; labels=transpose(Dmax), xlabel="t", ylabel=ob.name);
        end
        savefig(plt, string(savedir,unid,"/","convplot_",ob.name,"_",unid,".pdf"));
    end
end

function open_log(dt, tmax, convparams, method, machine, savedir, unid, name, params, obs, convobs, convcheck, kwargs...)
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
            writeprintln([f,f0], "\t method : $(method)")
            writeprintln([f,f0], "\t dt = $dt")
            writeprintln([f,f0], "\t tmax = $tmax")
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
            writeprintln([f,f0], "\t convparams : $convparams")
            #writeprint([f,f0], "\t options : ")
            #for (key, value) in kwargs
            #    writeprint([f,f0], String(key), " = ", value, ", ")
            #end
            #writeprintln([f,f0])
            writeprintln(f)
        end
    end
end

function error_log(savedir, unid)
    logname = string(savedir,"log-",gethostname(),".txt")
    open(logname, append=true) do f
        write(f, "[$(now())] => RUN <$unid> ERROR\n")
        write(f, "\t see $(unid)/$(unid).e for details\n")
        write(f, "\n")
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
        write(f, "\t total run time : $telapsed\n")
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

function save_plot(savedir, convcheck, unid, times, convdatadict, convparams, convobs)
    default(size = (800,600), reuse = true, title = unid, legend = true)
    numconv = length(convparams)
    labels = convcheck ? reshape(convparams, 1, numconv) : convparams
    for ob in filter(x->ndims(x)==0, convobs)
        if eltype(convdatadict[ob.name]) <: Complex
            plt = plot(convdatadict[ob.name]; labels=labels, xlabel="Re($(ob.name))", ylabel="Im($(ob.name))");
        else
            plt = plot(times, convdatadict[ob.name]; labels=labels, xlabel="t", ylabel=ob.name);
        end
        savefig(plt, string(savedir,unid,"/","convplot_",ob.name,"_",unid,".pdf"));
    end
end

macro LogParams(vars...)
    global names = string.(vars)
    quote
        len = length($(esc(vars)))
        vals = [$(esc.(vars)...)]
        [[names[i],vals[i]] for i=1:len]
    end
end

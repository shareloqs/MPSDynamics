#Deprecated#
function getchaincoeffs(nummodes, α, s, beta, ωc=1)
    matlabdir = ENV["MATDIR"]
    astr = paramstring(α, 2)
    bstr = paramstring(beta, 3)
    datfname = string("jacerg_","a",astr,"s$s","beta",bstr,".csv")
    chaincoeffs = readdlm(string(matlabdir,datfname),',',Float64,'\n')
    es = ωc .* chaincoeffs[:,1]
    ts = ωc .* chaincoeffs[1:end-1,2]
    c0 = ωc * sqrt(chaincoeffs[end,2]/pi)
    Nmax = length(es)
    if Nmax < nummodes
        throw(ErrorException("no data for nummodes > $Nmax"))
    end
    return [es[1:nummodes], ts[1:nummodes-1], c0]
end
##


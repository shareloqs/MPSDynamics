fdir = "/home/angus/Documents/Julia/mpos/"

α = 0.5
s = 1

coeffs = load(string(fdir,"coeffs_R1_s1_wc1.jld"))["coeffs"]
xmax, nmax = size(coeffs)
coeffs = [sqrt(2*α*(2*(n-1) + s + 1))*coeffs[x, n] for x=1:xmax, n=1:nmax]
chainparams = chaincoeffs_ohmic(500, α, s)


function electronphononmpo(Ne::Int, Nph::Int, d::Int, U=1.0, e=1.0, es=1.0, ts=1.0;
                           phonons = true,
                           couplingcoeffs = coeffs,
                           chainparams = chainparams
                           )

    N = Ne + Nph
    mpo = AutoMPO()
    for i=1:Ne
        mpo += (es, "Ntot", i)
    end
    for i=1:Ne-1
        mpo += (ts, "Cdagup", i, "Cup", i+1)
        mpo += (ts, "Cdagup", i+1, "Cup", i)
        mpo += (ts, "Cdagdn", i, "Cdn", i+1)
        mpo += (ts, "Cdagdn", i+1, "Cdn", i)
    end
    for i=1:Ne
        for j=i+1:Ne
            mpo += (e*(i-j)^(-2),"Ntot", i, "Ntot", j)
            #mpo += (10*exp(-norm(i-j)),"Ntot", i, "Ntot", j)
        end
    end
    for i=1:Ne
        mpo += (U, "Nup", i, "Ndn", i)
    end

    electronsites = [Index(4,"Electron") for i=1:Ne]

    if phonons
        ts = chainparams[2]
        es = chainparams[1]
        for i=Ne+1:N-1
            mpo += (ts[i-Ne], "A", i, "Adag", i+1)
            mpo += (ts[i-Ne], "Adag", i, "A", i+1)
        end
        for i=Ne+1:N
            mpo += (es[i-Ne], "N", i)
        end
        for i=1:Ne
            for j=Ne+1:N
                y = coeffs[i,j-Ne]
                mpo += (y, "Ntot", i, "A", j)
                mpo += (conj(y), "Ntot", i, "Adag", j)
            end
        end
        phononsites = [Index(d,"Boson") for i=1:Nph]
        sites = [electronsites..., phononsites...]
        return  MPOtoVector(MPO(mpo, sites))
    else
        return  MPOtoVector(MPO(mpo, electronsites))
    end
end
electronphononmpo(Ne::Int, U=1.0, e=1.0, es=1.0, ts=1.0) =
    electronphononmpo(Ne, 0, 1, U, e, es, ts; phonons = false)


#jldopen(string(fdir,"mpo_a0",".jld"), "w") do file
#    write(file, "H", H)
#end

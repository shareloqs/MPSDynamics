using MPSDynamics, Plots, LaTeXStrings, QuadGK, Revise

import MPSDynamics: disp, measuremodes, eigenchain, mpsembed!
import MPSDynamics: setleftbond, setrightbond, mpsrightnorm! 
import MPSDynamics: interleaved_tightbinding_mpo, productstatemps, tightbinding_mpo

#----------------------------
# Physical parameters
#----------------------------

N = 20  # number of chain sites
β = 2.0 # inverse temperature
μ = 0. # chemical potential
Ed = 0.6 
ϵd = Ed - μ # energy of the impurity

dt = 0.5
T = 30.0

dir = "/Users/ariva/Documents/fermions/"


# define the dispersion relation

function ϵ(x)
    return x
end


chainparams1 = chaincoeffs_fermionic(20, 2.0, 1.0; ϵ, save=true)
chainparams2 = chaincoeffs_fermionic(20, 2.0, 2.0; ϵ, save=true)

curdir = @__DIR__
dir_chaincoeff = abspath(joinpath(curdir, "../ChainOhmT/fermionicT/"))
chainparams2  = readchaincoeffs("$dir_chaincoeff/chaincoeffs.h5", β, 2.0) #filled
chainparams1  = readchaincoeffs("$dir_chaincoeff/chaincoeffs.h5", β, 1.0) #empty

#chainparams2 = readchaincoeffs(dir*"fermionic.h5", β, 2.0) 
#chainparams1 = readchaincoeffs(dir*"fermionic.h5", β, 1.0)

c1 = chainparams1[3]
c2 = chainparams2[3]

H = interleaved_tightbinding_mpo(N, ϵd, chainparams1, chainparams2)

ψ =  unitcol(2,2) # (0,1) filled impurity state
Tot = Any[]
push!(Tot, ψ)
for i in 1:(N)
    push!(Tot, unitcol(2,2))
    push!(Tot, unitcol(1,2))
end

A = productstatemps(physdims(H), state=Tot) # MPS
mpsembed!(A, 2)

dt = 0.5
T = 30.0

# Observables for the interleaved chain
ob1 = OneSiteObservable("system_occup", numb(2), 1)
ob2 = OneSiteObservable("folded_chain_occup", numb(2), (2,2N))

A, dat = runsim(dt, T, A, H;
                name = "Anderson impurity problem (folded chain)",
                #method = :TDVP2,
                #method = :TDVP1,
                method = :DTDVP,
                obs = [ob1, ob2], # for the interleaved chain
                #obs = [ob1, ob2, ob3], # for the double chain
                convobs = [ob1],
                params = @LogParams(N, ϵd, β, c1, c2),
                #convparams = [4],
                convparams = [0.1], # for the adaptive TDVP we set the precision that we want to obtain
                Dlim = 10,          # max bond dimension
                savebonddims = true, # we want to save the bond dimension
                verbose = false,
                save = true,
                savedir = dir,
                plot = true,
                );
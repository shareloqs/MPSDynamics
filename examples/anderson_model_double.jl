using MPSDynamics, Plots, LaTeXStrings, QuadGK

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


chainparams1 = chaincoeffs_fermionic(20, 2.0, 1.0; ϵ, save=false) # empty
chainparams2 = chaincoeffs_fermionic(20, 2.0, 2.0; ϵ, save=false) # filled

#=
curdir = @__DIR__
dir_chaincoeff = abspath(joinpath(curdir, "../ChainOhmT/fermionicT/"))
chainparams2  = readchaincoeffs("$dir_chaincoeff/chaincoeffs.h5", N, β, 2.0) #filled
chainparams1  = readchaincoeffs("$dir_chaincoeff/chaincoeffs.h5", N, β, 1.0) #empty
=#

c1 = chainparams1[3]
c2 = chainparams2[3]

H = tightbinding_mpo(N, ϵd, chainparams1, chainparams2)

ψ =  unitcol(2,2) # (0,1) filled impurity state
A = productstatemps(physdims(H), state=[fill(unitcol(2,2), N)..., ψ, fill(unitcol(1,2), N)...]) # MPS
mpsembed!(A, 4)

dt = 0.5
T = 30.0


# For the double chain  
ob1 = OneSiteObservable("chain1_filled_occup", numb(2), (1,N))
ob2 = OneSiteObservable("chain2_empty_occup", numb(2), (N+2, 2N+1))
ob3 = OneSiteObservable("system_occup", numb(2), N+1)

A, dat = runsim(dt, T, A, H;
                name = "Anderson impurity problem (folded chain)",
                #method = :TDVP2,
                #method = :TDVP1,
                method = :DTDVP,
                obs = [ob1, ob2, ob3], # for the double chain
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


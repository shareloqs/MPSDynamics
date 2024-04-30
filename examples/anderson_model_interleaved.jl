using MPSDynamics, Plots, LaTeXStrings, QuadGK, Revise
import MPSDynamics: disp, measuremodes, eigenchain, mpsembed!
import MPSDynamics: setleftbond, setrightbond, mpsrightnorm! 
import MPSDynamics: interleaved_tightbinding_mpo, productstatemps, tightbinding_mpo

#----------------------------
# Physical parameters
#----------------------------

N = 20      # number of chain sites
β = 2.0     # inverse temperature
μ = 0.      # chemical potential
Ed = 0.6    # energy of the impurity
ϵd = Ed - μ # energy of the impurity minus the chemical potential

#-----------------------------------------
# Dispersion relation and chain parameters
#-----------------------------------------

function ϵ(x)
    return x
end

chainparams1 = chaincoeffs_fermionic(20, 2.0, 1.0; ϵ, save=false) # empty
chainparams2 = chaincoeffs_fermionic(20, 2.0, 2.0; ϵ, save=false) # filled

#=
curdir = @__DIR__
dir_chaincoeff = abspath(joinpath(curdir, "../ChainOhmT/fermionicT/"))
chainparams2  = readchaincoeffs("$dir_chaincoeff/chaincoeffs.h5", β, 2.0) #filled
chainparams1  = readchaincoeffs("$dir_chaincoeff/chaincoeffs.h5", β, 1.0) #empty
=#

c1 = chainparams1[3]
c2 = chainparams2[3]

#-----------------------
# Simulation parameters
#-----------------------

dt = 0.5            # time step
T = 30.0            # simulation time
method = :DTDVP     # time-evolution method
Dmax = 2            # MPS max bond dimension
prec = 0.1          # precision for the adaptive TDVP

dir = "/Users/ariva/Documents/fermions/"

#---------------------------
# MPO and initial state MPS
#---------------------------

H = interleaved_tightbinding_mpo(N, ϵd, chainparams1, chainparams2)

ψ =  unitcol(2,2) # (0,1) filled impurity state
Tot = Any[]
push!(Tot, ψ)
for i in 1:(N)
    push!(Tot, unitcol(2,2))
    push!(Tot, unitcol(1,2))
end

A = productstatemps(physdims(H), state=Tot) # MPS
mpsembed!(A, 2) # embed the MPS in a manifold of bond dimension 2


#----------------------------------------------------
# Definition of observables for the interleaved chain
#----------------------------------------------------

ob1 = OneSiteObservable("system_occup", numb(2), 1)
ob2 = OneSiteObservable("folded_chain_occup", numb(2), (2,2N+1))

#-------------
# Simulation
#------------

A, dat = runsim(dt, T, A, H;
                name = "Anderson impurity problem (folded chain)",
                method = method,
                obs = [ob1, ob2], 
                convobs = [ob1],
                params = @LogParams(N, ϵd, β, c1, c2),
                convparams = [prec],   
                Dlim = Dmax,          
                savebonddims = true,   # we want to save the bond dimension
                verbose = false,
                save = true,
                savedir = dir,
                plot = true,
                );


#----------------------------
# Post-processing of the data
#----------------------------

occ_unfold = similar(dat["data/folded_chain_occup"])

# Loop through each element, modifying as necessary
for i in 1:size(dat["data/folded_chain_occup"], 1)
    for j in 1:size(dat["data/folded_chain_occup"], 2)
        if j % 2 != 0  # Check if the column index is odd
            occ_unfold[i, j] = 1 - dat["data/folded_chain_occup"][i, j]
        else
            occ_unfold[i, j] = dat["data/folded_chain_occup"][i, j]
        end
    end
end

#-------------
# Plots
#-------------


using MPSDynamics, Plots, LaTeXStrings, QuadGK

import MPSDynamics: disp, measuremodes, eigenchain, mpsembed!
import MPSDynamics: setleftbond, setrightbond, mpsrightnorm! 
import MPSDynamics: interleaved_tightbinding_mpo, productstatemps, tightbinding_mpo

#----------------------------
# Physical parameters
#----------------------------

N = 40      # number of chain sites
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

chainparams1 = chaincoeffs_fermionic(N, β, 1.0; ϵ, save=false) # empty
chainparams2 = chaincoeffs_fermionic(N, β, 2.0; ϵ, save=false) # filled

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
T = 10.0            # simulation time
method = :DTDVP     # time-evolution method
Dmax = 100           # MPS max bond dimension
prec = 0.0001       # precision for the adaptive TDVP

dir = "/Users/ariva/Documents/fermions/"

#---------------------------
# MPO and initial state MPS
#---------------------------

H = tightbinding_mpo(N, ϵd, chainparams1, chainparams2)

ψ =  unitcol(2,2) # (0,1) filled impurity state
A = productstatemps(physdims(H), state=[fill(unitcol(2,2), N)..., ψ, fill(unitcol(1,2), N)...]) # MPS
mpsembed!(A, 4)

#----------------------------------------------------
# Definition of observables for the double chain
#----------------------------------------------------

# For the double chain  
ob1 = OneSiteObservable("chain1_filled_occup", numb(2), (1,N))
ob2 = OneSiteObservable("chain2_empty_occup", numb(2), (N+2, 2N+1))
ob3 = OneSiteObservable("system_occup", numb(2), N+1)

#-------------
# Simulation
#------------

A, dat = runsim(dt, T, A, H;
                name = "Anderson impurity problem (folded chain)",
                method = method,
                obs = [ob1, ob2, ob3], 
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

#----------------
# Post-processing
#----------------

# Reshaping the vector to a column matrix and horizontal concatenation
system_occup_col = reshape(dat["data/system_occup"], :, 1)
occ = hcat(dat["data/chain1_filled_occup"]', system_occup_col)
occ = vcat(occ', dat["data/chain2_empty_occup"])

#-------------
# Plots
#-------------

# Plot the occupation of the chain sites
heatmap(
    collect(1:2*N+1),
    dat["data/times"],
    transpose(occ),  # Use the matrix form
    cmap = :coolwarm,
    aspect_ratio = :auto,
    xlabel = L"$N_{i,j}$ chain sites",
    ylabel = L"$t$",
    title = "Chain Occupation",
    colorbar = true,
    size = (700, 500)
)

# Plot the bond dimensions
heatmap(
    collect(1:2*N+2),
    dat["data/times"],
    transpose(dat["data/bonddims"]),
    cmap = :magma,
    aspect_ratio = :auto,
    xlabel = L"$N_{i,j}$ chain sites",
    ylabel = L"$t$",
    title = "Bond Dimensions",
    colorbar = true,
    size = (700, 500)
)

# Define indices for columns to be plotted
columns_to_plot = [1, 5, 10, 15, 20]

# Plot vertical slices for occupancy
p1 = plot(title = "Chain occupation")
for col in columns_to_plot
    plot!(p1, unfolded_occ_matrix[:, col], label = L"$t =$"*"$col", xlabel = L"$N_{i,j}$ chain sites", ylabel = "chain occupation")
end

# Plot vertical slices for bond dimensions
p2 = plot(title = "Bond Dimensions")
for col in columns_to_plot
    plot!(p2, unfolded_bonds_matrix[:, col], label = L"$t =$"*"$col", xlabel = L"$N_{i,j}$ chain sites", ylabel = L"$\chi$")
end

# Display the plots
plot(p1, p2, layout = (2, 1), size = (600, 800))

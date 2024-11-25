using MPSDynamics, Plots, LaTeXStrings
import MPSDynamics: productstatemps, tightbinding_mpo, mpsembed!

#----------------------------
# Physical parameters
#----------------------------

N = 40      # number of chain sites
β = 10.0     # inverse temperature
μ = 0.      # chemical potential
Ed = 0.3    # energy of the impurity
ϵd = Ed - μ # energy of the impurity minus the chemical potential

#-----------------------------------------
# Dispersion relation and chain parameters
#-----------------------------------------

function ϵ(x)
    return x
end

function J(x)
    return sqrt(1 - x^2) # semi-circular density of states
end

chainparams1 = chaincoeffs_fermionic(N, β, 1.0; ϵ, J, save=false) # empty
chainparams2 = chaincoeffs_fermionic(N, β, 2.0; ϵ, J, save=false) # filled

c1 = chainparams1[3]
c2 = chainparams2[3]  

#=
curdir = @__DIR__
dir_chaincoeff = abspath(joinpath(curdir, "../ChainOhmT/fermionicT/"))
chainparams2  = readchaincoeffs("$dir_chaincoeff/chaincoeffs.h5", β, 2.0) #filled
chainparams1  = readchaincoeffs("$dir_chaincoeff/chaincoeffs.h5", β, 1.0) #empty
=#

#-----------------------
# Simulation parameters
#-----------------------

dt = 0.25            # time step
T = 15.0            # simulation time
method = :DTDVP     # time-evolution method
Dmax = 150           # MPS max bond dimension
prec = 0.00001       # precision for the adaptive TDVP

#---------------------------
# MPO and initial state MPS
#---------------------------

H = tightbinding_mpo(N, ϵd, chainparams1, chainparams2)

ψ =  unitcol(2,2) # (0,1) filled impurity state
A = productstatemps(physdims(H), state=[fill(unitcol(2,2), N)..., ψ, fill(unitcol(1,2), N)...]) # MPS
mpsembed!(A, 2)

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
                save = false,
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

# Plot the system occupation    
p1 = plot(
    dat["data/times"],
    dat["data/system_occup"],
    xlabel = L"$t$",
    ylabel = L"$n_d$",
    title = "System Occupation",
    size = (700, 500)
)

# Plot the occupation of the chain sites
p2 = heatmap(
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
p3 = heatmap(
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
p4 = plot(title = "Chain occupation")
for col in columns_to_plot
    plot!(p4, occ[:, col], label = L"$t =$"*"$col", xlabel = L"$N_{i,j}$ chain sites", ylabel = "chain occupation")
end

# Plot vertical slices for bond dimensions
p5 = plot(title = "Bond Dimensions")
for col in columns_to_plot
    plot!(p5, dat["data/bonddims"][:, col], label = L"$t =$"*"$col", xlabel = L"$N_{i,j}$ chain sites", ylabel = L"$\chi$")
end

# Display the plots
plot(p2, p3, p4, p5, p1, layout = (3, 2), size = (1400, 1200))

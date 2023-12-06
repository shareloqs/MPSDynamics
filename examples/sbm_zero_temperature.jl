"""
    Example of a zero-temperature Spin-Boson Model with an hard cut-off Ohmic spectral density J(ω) = 2αω when ω < ωc and 0 otherwise

    The dynamics is simulated using the T-TEDOPA method that maps the normal modes environment into a non-uniform tight-binding chain.

    H = \frac{ω_0}{2} σ_z + Δ σ_x + c_0 σ_x(b_0^\dagger + b_0) + \sum_{i=0}^{N-1} t_i (b_{i+1}^\dagger b_i +h.c.) + \sum_{i=0}^{N-1} ϵ_i b_i^\dagger b_i  
"""

using MPSDynamics, Plots, LaTeXStrings

#----------------------------
# Physical parameters
#----------------------------

d = 6 # number of Fock states of the chain modes

N = 30 # length of the chain

α = 0.1 # coupling strength

Δ = 0.0 # tunneling 

ω0 = 0.2 # TLS gap

s = 1 # ohmicity

cpars = chaincoeffs_ohmic(N, α, s) # chain parameters, i.e. on-site energies ϵ_i, hopping energies t_i, and system-chain coupling c_0

#-----------------------
# Simulation parameters
#-----------------------

dt = 0.5 # time step

tfinal = 30.0 # simulation time

method = :TDVP1 # time-evolution method

D = [2,4,6] # MPS bond dimension

#---------------------------
# MPO and initial state MPS
#---------------------------

H = spinbosonmpo(ω0, Δ, d, N, cpars) # MPO representation of the Hamiltonian

ψ = unitcol(1,2) # Initial up-z system state 

A = productstatemps(physdims(H), state=[ψ, fill(unitcol(1,d), N)...]) # MPS representation of |ψ>|Vacuum>

#---------------------------
# Definition of observables
#---------------------------

ob1 = OneSiteObservable("sz", sz, 1)

ob2 = OneSiteObservable("chain mode occupation", numb(d), (2,N+1))

ob3 = TwoSiteObservable("SXdisp", sx, disp(d), [1], collect(2:N+1))

#-------------
# Simulation
#------------

A, dat = runsim(dt, tfinal, A, H;
                name = "ohmic spin boson model",
                method = method,
                obs = [ob2,ob3],
                convobs = [ob1],
                params = @LogParams(N, d, α, Δ, ω0, s),
                convparams = D,
                verbose = false,
                save = true,
                plot = true,
                );

#----------
# Plots
#----------

plot(dat["data/times"], dat["convdata/sz"], label=["Dmax = 2" "Dmax = 4" "Dmax = 6"], xlabel=L"t",ylabel=L"\sigma_z")

heatmap(dat["data/times"], collect(1:N), abs.(dat["data/SXdisp"][1,:,:]), xlabel=L"t",ylabel="chain mode")

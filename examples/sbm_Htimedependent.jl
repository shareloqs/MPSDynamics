#=
    Example of a time-dependent Hamitloninan on a zero-temperature Spin-Boson Model with an hard cut-off Ohmic spectral density J(ω) = 2αω when ω < ωc and 0 otherwise

    The dynamics is simulated using the T-TEDOPA method that maps the normal modes environment into a non-uniform tight-binding chain.

    H = \\frac{ω_0}{2} σ_z + Δ σ_x + σ_x ϵ sin(ωdrive t) + c_0 σ_x(b_0^\\dagger + b_0) + \\sum_{i=0}^{N-1} t_i (b_{i+1}^\\dagger b_i +h.c.) + \\sum_{i=0}^{N-1} ϵ_i b_i^\\dagger b_i  
=#

using MPSDynamics, Plots, LaTeXStrings, LinearAlgebra

import MPSDynamics: disp

#----------------------------
# Physical parameters
#----------------------------

d = 4 # number of Fock states of the chain modes

N = 60 # length of the chain

α = 0.005 # coupling strength

Δ = 0.0 # tunneling 

ω0 = 0.8 # TLS gap

s = 1 # ohmicity

cpars = chaincoeffs_ohmic(N, α, s) # chain parameters, i.e. on-site energies ϵ_i, hopping energies t_i, and system-chain coupling c_0

Trabi = 30.0 # Rabi period of the drive

ϵ = 2*pi / Trabi # Intensity of the drive

ωdrive = ω0 # Frequency of the drive

#-----------------------
# Simulation parameters
#-----------------------

dt = 0.5 # time step

tfinal = 100.0 # simulation time

method = :TDVP1 # time-evolution method

D = [6] # MPS bond dimension

#---------------------------
# MPO time-dependent
#---------------------------

timelist = collect(0:dt:tfinal)
numsteps = length(timelist)-1

#### #=
    #Example of a Ht = σ_x ϵ sin(ωdrive t) drive on the TLS
    #To comment if the example of (a+a^\dagger) drive is used

Ndrive = 1 #Number of the site on which the drive is applied. Here Ndrive=1 for the TLS
Ht = [ϵ*sx*sin(ωdrive*tstep) for tstep in timelist] # Time-dependent Hamiltonian term

#### =#

#=
   # Example of a Ht = (a_{Ndrive_star}+a_{Ndrive_star}^\dagger) ϵ sin(ωdrive t) drive on the (Ndrive_star)th mode. 
   # The operators have to be transformed to drive the chain modes instead of the mode of the initial Hamiltonian
   # This example assumes the same dimensions for all Ndrive and that applies to a chain at the right of the system.
   # If it is not the case, the run_1TDVP loop involved when timedep=true has to be modified
   # To comment if the example of σ_x drive is used

Ndrive_star = 10
Ndrive = collect(2:N+1)

# Construct the chain Hamiltonian to write the star - chain transition matrix 
t_chain=cpars[2][1:N-1] ; e_chain=cpars[1][1:N]
hmat_chain = MPSDynamics.diagm(0=>e_chain, 1=>t_chain, -1=>t_chain)
U_chain  = eigen(hmat_chain).vectors

print("\n Driving of the mode of frequency : ",MPSDynamics.eigenchain(cpars, nummodes=N).values[Ndrive_star])

Ht = Vector{Any}(undef,N+1)

# Construct the chain operators that correspond to star operators for (a_{Ndrive_star}+a_{Ndrive_star}^\dagger
for i in Ndrive
    local anih_Ndrive = (U_chain[i-1,Ndrive_star]*anih(d))
    local crea_Ndrive = (conj(U_chain[i-1,Ndrive_star])*crea(d))
    local Ht_Ndrive = [ϵ*(anih_Ndrive+crea_Ndrive)*sin(ωdrive*tstep) for tstep in timelist]
    Ht[i] = Ht_Ndrive
end

=#
#---------------------------
# MPO and initial state MPS
#---------------------------

H = spinbosonmpo(ω0, Δ, d, N, cpars) # MPO representation of the Hamiltonian

ψ = unitcol(2,2) # Initial down-z system state 

A = productstatemps(physdims(H), state=[ψ, fill(unitcol(1,d), N)...]) # MPS representation of |ψ>|Vacuum>

#---------------------------
# Definition of observables
#---------------------------

ob1 = OneSiteObservable("sz", sz, 1)

ob2 = TwoSiteObservable("cdagc", crea(d), anih(d), collect(2:N+1), collect(2:N+1))

ob3 = TwoSiteObservable("SXdisp", sx, disp(d), [1], collect(2:N+1))


#-------------
# Simulation
#------------

A, dat = runsim(dt, tfinal, A, H;
                name = "Driving field on ohmic spin boson model",
                method = method,
                obs = [ob1, ob2, ob3],
                convobs = [ob1],
                params = @LogParams(N, d, α, Δ, ω0, s),
                convparams = D,
                timedep = true, # the Hamiltonian is time dependent
                Ndrive = Ndrive, # site(s) of the MPS/MPO that are driven
                Htime = Ht, # list of time-dependent terms
                verbose = false,
                save = true,
                plot = true,
                );

#----------
# Plots
#----------

plot(dat["data/times"], dat["data/sz"], label="Dmax = $(D...)", xlabel=L"t",ylabel=L"\sigma_z", title="")

bath_occup = mapslices(X -> MPSDynamics.measuremodes(X, cpars[1], cpars[2]), dat["data/cdagc"], dims = [1,2])
omeg = MPSDynamics.eigenchain(cpars, nummodes=N).values

plot(omeg, bath_occup[:, :, end], lw=4, xlabel=L"\omega", ylabel=L"\langle n^b_\omega \rangle",
title="Mode occupation in the extended bath at final time")
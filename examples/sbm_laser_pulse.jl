#=
    Example of a time-dependent Hamitloninan on a zero-temperature Spin-Boson Model with an hard cut-off Ohmic spectral density J(ω) = 2αω when ω < ωc and 0 otherwise

    The dynamics is simulated using the T-TEDOPA method that maps the normal modes environment into a non-uniform tight-binding chain.

    H = \\frac{ω_0}{2} σ_z + Δ σ_x + σ_x ϵ sin(ωdrive t) + c_0 σ_x(b_0^\\dagger + b_0) + \\sum_{i=0}^{N-1} t_i (b_{i+1}^\\dagger b_i +h.c.) + \\sum_{i=0}^{N-1} ϵ_i b_i^\\dagger b_i  
=#

using MPSDynamics, Plots, LaTeXStrings

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

Ndrive = 1 #Number of the site on which the drive is applied
#-----------------------
# Simulation parameters
#-----------------------

dt = 0.5 # time step

tfinal = 100.0 # simulation time

method = :TDVP1 # time-evolution method

D = [6] # MPS bond dimension

#---------------------------
# MPO and initial state MPS
#---------------------------

numsteps = length(collect(0:dt:tfinal))-1
timelist = [(i-1)*dt for i=1:numsteps+1]

Ht = [ϵ*sx*sin(ωdrive*tstep) for tstep in timelist] # Time dependent Hamiltonian term MPO

H = spinbosonmpo(ω0, Δ, d, N, cpars) # MPO representation of the Hamiltonian

ψ = unitcol(2,2) # Initial low-z system state 

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
                name = "Driving field on ohmic spin boson model",
                method = method,
                obs = [ob1],
                convobs = [ob1],
                params = @LogParams(N, d, α, Δ, ω0, s),
                convparams = D,
                timedep = true,
                Ndrive = Ndrive,
                Htime = Ht,
                verbose = false,
                save = true,
                plot = true,
                );

#----------
# Plots
#----------

plot(dat["data/times"], dat["data/sz"], label=["Dmax = $D"], xlabel=L"t",ylabel=L"\sigma_z", title="")
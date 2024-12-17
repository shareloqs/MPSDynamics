#=
    Example of a correlated environment model at finite or zero temperature with a hard cut-off Ohmic spectral density J(ω) = 2α(ω^s)/(ω_c^(s-1)) when 0 <= ω <= ωc and 0 otherwise

    A correlated environment is a global environment interacting with a multi-site system.

    H = \sum_n \frac{E_n}{2}σ_n^z + J_n (σ_n^- σ_{n+1}^+ + h.c.) + \int_\mathbb{R} \omega_k \a^\dagger_k \a_k dk + \sum_n 0.5*(1 + σ_n^z)\int_\mathbb{R}( \sqrt{J(\omega_k)} \exp(i k r_n) a_k + h.c.) dk
=#

using MPSDynamics, Plots, LaTeXStrings

import MPSDynamics: correlatedenvironmentmpo, mpsembed!

const ∞  = Inf

#----------------------------
# Physical parameters
#----------------------------

N = 2 # number of system sites

R = [i*10 for i=0:N-1] # positions of the sites

E = [0., 0.] # sites' energies

J = [0.25] # sites' tunneling

P = [1 0;0 0] # local excited state projector

As = [P for i=1:N] # system operators of the interaction Hamiltonian H_int = \sum_s A_s \otimes B_s

d = 5 # number of Fock states of the chain modes

M = 15 # length of the chain

α = 0.1 # coupling strength

s = 1 # ohmicity

ωc = 1 # Cut-off of the spectral density J(ω)

c = 1 # phonon speed of sound

#β = 10 # Thermalized environment

β = ∞ # Case zero temperature T=0, β → ∞

fnamecc = "./couplings.csv" # path to the file that (will) contains the system/environment coupling constants 

#----------------------------
# Ohmic spectral density
#----------------------------

if β == ∞
    cpars = chaincoeffs_ohmic(M, α, s; ωc=ωc)  # chain parameters, i.e. on-site energies ϵ_i, hopping energies t_i, and system-chain coupling c_0
else
    cpars = chaincoeffs_finiteT(M, β; α=α, s=s, J=nothing, ωc=ωc, mc=4, mp=0, AB=nothing, iq=1, idelta=2, procedure=:Lanczos, Mmax=5000, save=true)  # chain parameters, i.e. on-site energies ϵ_i, hopping energies t_i, and system-chain coupling c_0
    #=#If cpars is stored in "../ChainOhmT/ohmicT" 
    curdir = @__DIR__
    dir_chaincoeff = abspath(joinpath(curdir, "../ChainOhmT/ohmicT"))
    cpars  = readchaincoeffs("$dir_chaincoeff/chaincoeffs.h5", N, α, s, β) # chain parameters, i.e. on-site energies ϵ_i, hopping energies t_i, and system-chain coupling c_0
    =#
end

#-----------------------
# Simulation parameters
#-----------------------

dt = 0.5 # time step

tfinal = 30.0 # simulation time

method = :DTDVP # time-evolution method

eps = [1e-3] # projection error threshol for DTDVP

#---------------------------
# MPO and initial state MPS
#---------------------------

Hs = MPSDynamics.multisitempo(N, E, J, As) # system Hamiltonian

Hc = correlatedenvironmentmpo(R, M, d; chainparams=cpars, fnamecc=fnamecc, s=s, α=α, ωc=ωc, c_phonon=c, β=β) # environment Hamiltonian

H = [Hs..., Hc...] # total Hamiltonian

# Initial electronic system in a superposition of 1 and 2

ψ = [unitcol(1,2), unitcol(2,2)] # system initial state with an excitation on the first site

A = productstatemps(physdims(H), state=[ψ..., fill(unitcol(1,d), 2M)...]) # MPS representation of |ψ>|Vacuum>

# In the single-excitation subspace, a superposition of an excitation on the first site and an excitation on the second site is an entangled state |qubit 1 up>|qubit 2 down> + |qubit 1 down>|qubit 2 up>
# We here construct this initial entangled state by hand

#mpsembed!(A,2)

#temp = zeros(ComplexF64, 1, 2, 2)
#temp[:,:,1] = [1/sqrt(2) 0]
#temp[:,:,2] = [0 1/sqrt(2)]
#A[1] = temp
#temp = zeros(ComplexF64, 2, 2, 2)
#temp[:,:,1] = [0 0;1 0]
#temp[:,:,2] = [1 0;0 0]
#A[2] = temp

#---------------------------
# Definition of observables
#---------------------------

ob1 = OneSiteObservable("P1", P, 1)

ob2 = OneSiteObservable("P2", P, 2)

ob3 = OneSiteObservable("chain mode occupation", numb(d), (N+1,2M+N))

ob4 = TwoSiteObservable("coherence", crea(2), anih(2), 1,2)

#-------------
# Simulation
#------------

A, dat = runsim(dt, tfinal, A, H;
                name = "correlated environment",
                method = method,
                obs = [ob1, ob2, ob3,ob4],
                convobs = [ob1],
                params = @LogParams(E, J, As, N, M, d, α, s, β),
                convparams = eps,
                verbose = false,
                save = true,
                plot = true,
                savebonddims = true
                );

#-------------
# Plots
#-------------

method == :DTDVP && heatmap(dat["data/times"], collect(0:N+2M), dat["data/bonddims"], xlabel=L"\omega_c t",ylabel="bond index", cbar_title=L"\mathrm{Bond dimension}\ \chi")

heatmap(collect(-M:M-1), dat["data/times"], vcat(dat["data/chain mode occupation"][M:-1:1,:], dat["data/chain mode occupation"][M+1:2M,:])', xlabel=L"\mathrm{chain mode,}\ n", ylabel=L"\omega_c t", cbar_title=L"\langle\hat{b}_n^\dagger\hat{b}_n\rangle(t)")

plot(dat["data/times"], [dat["data/P1"] dat["data/P2"]], xlabel=L"\omega_c t", ylabel="Population", label=["Site 1" "Site 2"], linewidth=4, titlefontsize=16, legend=:best, legendfontsize=16, xguidefontsize=16, yguidefontsize=16, tickfontsize=10)

#pup = 0.5 .+ real.(dat["data/coherence"]) # Upper level population in the degenerate case
## For the general case one can either solve a second degree equation involving the population of one of the sites and the coherence and pick at each time the physical solution, or one can compute the density matrix and diagonalize it at each time

#plot(dat["data/times"], pup, xlabel=L"\omega_c t", ylabel="Population", label="Upper Level", linewidth=4, titlefontsize=16, legend=:best, legendfontsize=16, xguidefontsize=16, yguidefontsize=16, tickfontsize=10)

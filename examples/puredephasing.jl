#=
    Example of a Pure Dephasing Model at finite or zero temperature with an hard cut-off Ohmic spectral density J(ω) = 2α(ω^s)/(ω_c^(s-1)) when ω < ωc and 0 otherwise#

    The dynamics is simulated using the T-TEDOPA method that maps the normal modes environment into a non-uniform tight-binding chain.

    H = \frac{ω0}{2} σ_z +  \frac{σ_z}{2} c_0 (b_0^\dagger + b_0) + \sum_{i=0}^{N-1} t_i (b_{i+1}^\dagger b_i +h.c.) + \sum_{i=0}^{N-1} ϵ_i b_i^\dagger b_i  
=#

using MPSDynamics, Plots, LaTeXStrings, QuadGK

import MPSDynamics: chaincoeffs_ohmic, puredephasingmpo

const ∞  = Inf
#----------------------------
# Physical parameters
#----------------------------

d = 10 # number of Fock states of the chain modes

N = 30 # length of the chain

α = 0.01 # coupling strength

ω0 = 0.008 # TLS gap

s = 1 # ohmicity

ωc = 0.035 # Cut-off of the spectral density J(ω)

#β = 100 # Thermalized environment

β = ∞ # Case zero temperature T=0, β → ∞

#----------------------------
# Ohmic spectral density
#----------------------------

if β == ∞
    cpars = chaincoeffs_ohmic(N, α, s; ωc=ωc)  # chain parameters, i.e. on-site energies ϵ_i, hopping energies t_i, and system-chain coupling c_0
else
    cpars = chaincoeffs_finiteT(N, β; α=α, s=s, J=nothing, ωc=ωc, mc=4, mp=0, AB=nothing, iq=1, idelta=2, procedure=:Lanczos, Mmax=5000, save=true)  # chain parameters, i.e. on-site energies ϵ_i, hopping energies t_i, and system-chain coupling c_0
    #=#If cpars is stored in "../ChainOhmT/ohmicT" 
    curdir = @__DIR__
    dir_chaincoeff = abspath(joinpath(curdir, "../ChainOhmT/ohmicT"))
    cpars  = readchaincoeffs("$dir_chaincoeff/chaincoeffs.h5", N, α, s, β) # chain parameters, i.e. on-site energies ϵ_i, hopping energies t_i, and system-chain coupling c_0
    =#
end


#-----------------------
# Simulation parameters
#-----------------------

dt = 1.0 # time step

tfinal = 10.0 # simulation time

method = :TDVP1 # time-evolution method

#method = :DTDVP # time-evolution method

D = 2 # MPS bond dimension

#---------------------------
# MPO and initial state MPS
#---------------------------

H = puredephasingmpo(ω0, d, N, cpars)

# Initial electronic system in a superposition of 1 and 2
ψ = zeros(2)
ψ[1] = 1/sqrt(2)
ψ[2] = 1/sqrt(2)


A = productstatemps(physdims(H), state=[ψ, fill(unitcol(1,d), N)...]) # MPS representation of |ψ>|Vacuum>

#---------------------------
# Definition of observables
#---------------------------

ob1 = OneSiteObservable("sz", sz, 1)


#-------------
# Simulation
#------------

A, dat = runsim(dt, tfinal, A, H, prec=1E-4;
                name = "pure dephasing model with temperature",
                method = method,
                obs = [ob1],
                convobs = [ob1],
                params = @LogParams(ω0, N, d, α, s),
                convparams = D,
                reduceddensity=true,
                verbose = false,
                save = true,
                plot = true,
                );



#----------
# Analytical results at specified temperature 
# (see The Theory of Open Quantum System, H.-P. Breuer & F. Petruccione 2002, Chapter 4)
#----------

Johmic(ω,s) = (2*α*ω^s)/(ωc^(s-1))

time_analytical = LinRange(0.0,tfinal,Int(tfinal))

Γohmic(t) = - quadgk(x -> Johmic(x,s)*(1 - cos(x*t))*coth(β*x/2)/x^2, 0, ωc)[1]

Decoherence_ohmic(t) = 0.5*exp(Γohmic(t))

#-------------
# Plots
#------------

ρ12 = abs.(dat["data/Reduced ρ"][1,2,:])

plot(time_analytical, t->Decoherence_ohmic(t), label="Analytics", title=L"Pure Dephasing, Ohmic $s=%$s$, $\beta = %$β ~\mathrm{K}$", linecolor=:black, xlabel="Time (arb. units)", ylabel=L"Coherence $|\rho_{12}(t)|$", linewidth=4, titlefontsize=16, legend=:best, legendfontsize=16, xguidefontsize=16, yguidefontsize=16, tickfontsize=10)

plot!(dat["data/times"], ρ12, lw=4, ls=:dash, label="Numerics")


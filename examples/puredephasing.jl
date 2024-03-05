#=
    Example of a Pure Dephasing Model at zero temperature and with an arbitrary temperature.  with an hard cut-off Ohmic spectral density J(ω) = 2αω when ω < ωc and 0 otherwise#

    The dynamics is simulated using the T-TEDOPA method that maps the normal modes environment into a non-uniform tight-binding chain.

    H = \frac{ΔE}{2} σ_z +  \frac{σ_z}{2} c_0 (b_0^\dagger + b_0) + \sum_{i=0}^{N-1} t_i (b_{i+1}^\dagger b_i +h.c.) + \sum_{i=0}^{N-1} ϵ_i b_i^\dagger b_i  
=#

using MPSDynamics, Plots, LaTeXStrings, QuadGK

#----------------------------
# Physical parameters
#----------------------------

ΔE = 0.008 # Energy of the electronic states

d = 5 # number of Fock states of the chain modes

N = 30 # length of the chain

α = 0.01 # coupling strength

s = 1 # ohmicity

ωc = 0.035 # Cut-off of the spectral density J(ω)

cpars = chaincoeffs_ohmic(N, α, s; ωc=ωc) # chain parameters, i.e. on-site energies ϵ_i, hopping energies t_i, and system-chain coupling c_0

#-----------------------
# Simulation parameters
#-----------------------

dt = 1.0 # time step

tfinal = 300.0 # simulation time

method = :TDVP1 # time-evolution method

D = 6 # MPS bond dimension

#---------------------------
# MPO and initial state MPS
#---------------------------

H = puredephasingmpo(ΔE, d, N, cpars)

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

A, dat = runsim(dt, tfinal, A, H;
                name = "pure dephasing model at zero temperature",
                method = method,
                obs = [ob1],
                convobs = [ob1],
                params = @LogParams(ΔE, N, d, α, s),
                convparams = D,
                reduceddensity=true,
                verbose = false,
                save = true,
                plot = true,
                );



#----------
# Analytical results at zero temperature 
# (see The Theory of Open Quantum System, H.-P. Breuer & F. Petruccione 2002, Chapter 4)
#----------

Johmic(ω,s) = (2*α*ω^s)/(ωc^(s-1))
f(ω,t) = (1 - cos(ω*t))/ω^2
time_analytical = LinRange(0.0,tfinal,Int(tfinal))

Γohmic(t) = - quadgk(x -> Johmic(x,s)*f(x,t), 0, ωc)[1]
Decoherence_ohmic(t) = 0.5*exp(Γohmic(t))
ListDecoherence_ohmic = [Decoherence_ohmic(t) for t in time_analytical]

#-------------
# Plots
#------------

ρ12 = sqrt.(real(dat["data/Reduced ρ"][1,2,:]).^2 .+ imag(dat["data/Reduced ρ"][1,2,:]).^2 )

(plot(time_analytical, ListDecoherence_ohmic1, label="Analytics", title=L"Pure Dephasing, Ohmic $s=%$s \, ,\, T=0K$", linecolor =:black, xlabel="Time (arb. units)",ylabel="Coherence Amplitude", linewidth=4, titlefontsize=16, legend=:best, legendfont=16, xguidefontsize=16, yguidefontsize=16, tickfontsize=10))#,ylims=(0.25,0.5)))#,xticks=(0:2.5:10))))
display(plot!(dat["data/times"],ρ12,lw=4,ls=:dash,label="Numerics"))
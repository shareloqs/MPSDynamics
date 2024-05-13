using MPSDynamics, Plots, LaTeXStrings, QuadGK, LinearAlgebra

import MPSDynamics: measuremodes, measurecorrs, mpsembed!, eigenchain

const ∞  = Inf
#----------------------------
# Physical parameters
#----------------------------

d = 10      # number of Fock states of the chain modes
N = 30      # length of the chain
α = 0.01    # coupling strength
ω0 = 0.008  # TLS gap
s = 1       # ohmicity
ωc = 0.035  # Cut-off of the spectral density J(ω)
#β = 100    # Thermalized environment
β = ∞       # Case zero temperature T=0, β → ∞

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

#method = :TDVP1         # time-evolution method
#conv = 2                # MPS bond dimension
#dt = 1.0                # time step
#tfinal = 10.0           # simulation time

method = :TDVP1         # time-evolution method
conv = 3                # bond dimension for the TDVP1
dt = 0.5                # time step
tfinal = 100.0           # simulation time
Tsteps = Int(tfinal / dt)

#---------------------------
# MPO and initial state MPS
#---------------------------

H = puredephasingmpo(ω0, d, N, cpars)

# Initial electronic system in a superposition of 1 and 2
ψ = zeros(2)
ψ[1] = 1/sqrt(2)
ψ[2] = 1/sqrt(2)


A = productstatemps(physdims(H), state=[ψ, fill(unitcol(1,d), N)...]) # MPS representation of |ψ>|Vacuum>
mpsembed!(A, 2)
#---------------------------
# Definition of observables
#---------------------------

ob1 = OneSiteObservable("sz", sz, 1)
ob2 = OneSiteObservable("sx", sx, 1)
ob3 = OneSiteObservable("chain_mode_occupation", numb(d), (2,N+1))
ob4 = OneSiteObservable("c", crea(d), collect(2:N+1))
ob5 = OneSiteObservable("cdag", crea(d), collect(2:N+1))

ob6 = TwoSiteObservable("cdagc", crea(d), anih(d), collect(2:N+1), collect(2:N+1))
ob7 = TwoSiteObservable("cdagcdag", crea(d), crea(d), collect(2:N+1), collect(2:N+1))
ob8 = TwoSiteObservable("cc", anih(d), anih(d), collect(2:N+1), collect(2:N+1))

#-------------
# Simulation
#------------

A, dat = runsim(dt, tfinal, A, H, prec=1E-4;
                name = "Bath observables in the pure dephasing model",
                method = method,
                obs = [ob1, ob2, ob3, ob4, ob5, ob6, ob7, ob8],
                convobs = [ob1],
                params = @LogParams(ω0, N, d, α, s, ψ),
                convparams = conv,
                Dlim = 100,
                reduceddensity = true,
                verbose = false,
                save = false,
                plot = true,
                );

#---------------------------
# Post-processing
#---------------------------

T = length(dat["data/times"])

constr = Array{ComplexF64}(undef, N, N, T)
destr = Array{ComplexF64}(undef, N, N, T)
for t in 1:T
    constr[:,:,t] = diagm(0 => dat["data/cdag"][:,t])
    destr[:,:,t] = diagm(0 => dat["data/c"][:,t])
end

omeg = eigenchain(cpars, nummodes=N).values
bath_occup = mapslices(X -> measuremodes(X, cpars[1], cpars[2]), dat["data/cdagc"], dims = [1,2])
cdag_average = mapslices(X -> measuremodes(X, cpars[1], cpars[2]), constr, dims = [1,2])
c_average = mapslices(X -> measuremodes(X, cpars[1], cpars[2]), destr, dims = [1,2])
cc_average = mapslices(X -> measurecorrs(X, cpars[1], cpars[2]), dat["data/cc"], dims = [1,2])
cdagcdag_average = mapslices(X -> measurecorrs(X, cpars[1], cpars[2]), dat["data/cdagcdag"], dims = [1,2])

correlations_c = [
    cc_average[i, j, t] - c_average[i, 1, t] .* c_average[j, 1, t]
    for i in 1:size(cc_average, 1), j in 1:size(cc_average, 2), t in 1:size(cc_average, 3)
]
correlations_cdag = [
    cdagcdag_average[i, j, t] - cdag_average[i, 1, t] .* cdag_average[j, 1, t]
    for i in 1:size(cdagcdag_average, 1), j in 1:size(cdagcdag_average, 2), t in 1:size(cdagcdag_average,3)
]
#--------------------
# Analytical results 
#--------------------

Johmic(ω,s) = (2*α*ω^s)/(ωc^(s-1))

time_analytical = LinRange(0.0, tfinal, Int(tfinal))

Γohmic(t) = - quadgk(x -> Johmic(x,s)*(1 - cos(x*t))*coth(β*x/2)/x^2, 0, ωc)[1]

Decoherence_ohmic(t) = 0.5 * exp(Γohmic(t))


α_theo = 0.25 * α
function Jtherm(x)
    if 1 >= x >= 0
        return +α_theo * abs(x)^s * (1 + coth(β*0.5*x))
    elseif -1 <= x <= 0
        return -α_theo * abs(x)^s * (1 + coth(β*0.5*x))
    else
        return 0
    end
end

bath_occup_analytical(ω, t) = abs(Jtherm(ω))/(ω^2)*2*(1-cos(ω*t)) 

#-------------
# Plots
#------------

ρ12 = abs.(dat["data/Reduced ρ"][1,2,:])

p1 = plot(time_analytical, t->Decoherence_ohmic(t), label="Analytics", title=L"Pure Dephasing, Ohmic $s=%$s$, $\beta = %$β ~\mathrm{K}$", linecolor=:black, xlabel="Time (arb. units)", ylabel=L"Coherence $|\rho_{12}(t)|$", linewidth=4, titlefontsize=16, legend=:best, legendfontsize=16, xguidefontsize=16, yguidefontsize=16, tickfontsize=10)
p1 = plot!(dat["data/times"], ρ12, lw=4, ls=:dash, label="Numerics")

cumul = [bath_occup_analytical(omeg[i], tfinal)*(omeg[i+1]-omeg[i]) for i in 1:(length(omeg)-1)]

p2 = plot(omeg[1:length(omeg)-1], cumul,
             xlabel=L"\omega", ylabel=L"\langle n^b_\omega \rangle", label="Analytics",
             title="Mode occupation")
p2 = plot!(omeg, bath_occup[:, :, Tsteps], label="Numerics")

p3 = heatmap(omeg, omeg, abs.(real.(correlations_cdag[:,:,Tsteps]) .+ im*imag.(correlations_cdag[:,:,Tsteps])), 
            xlabel=L"\omega",
            ylabel=L"\omega", title="Environmental correlations")

plot(p1, p2, p3, layout = (2, 2), size = (1400, 1200))

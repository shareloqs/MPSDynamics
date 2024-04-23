# Pure-Dephasing

## Context 

The Pure-Dephasing Model describes a two-level system interacting linearly with an environment characterised by a spectral density (SD) ``J(\omega)``. The coupling only acts on diagonal terms through the ``\sigma_z`` operator. The Hamiltonian reads

```math
        \hat{H} = \frac{\omega_0}{2}\hat{\sigma}_z + \int_0^{+\infty}\omega \hat{a}^\dagger_\omega\hat{a}_\omega \mathrm{d}\omega + \frac{\hat{\sigma}_z}{2}\int_0^{+\infty}\sqrt{J(\omega)}(\hat{a}_\omega + \hat{a}^\dagger_\omega)\mathrm{d}\omega
```

Although this interaction will not change the population of the two-level system, the coherences between the two states will vary due to the environment. Introducing the two-level system reduced density matrix ``\rho_{ij}(t)`` with ``i,j \in (0,1)``, the diagonal terms ``\rho_{ii}(t)`` are the populations of the states and the anti-diagonal terms ``\rho_{ij}(t)`` with ``i \neq j `` are the coherences between the two states. The effect of the ``\sigma_z`` bath interaction is to decouple the two states ``|0\rangle`` and ``|1\rangle``. 

The density matrix can be calculated within the MPS formalism with the MPO ``|\psi\rangle \langle \psi|``. Tracing out the environment leads to the reduced density matrix. This can be done with the function [`MPSDynamics.rhoreduced_1site`](@ref). 

An analytical formula can be found for the decoherence function ``\Gamma(t)``, taking into account the SD as well as the temperature of the environment [^breuer]

```math
    \Gamma(t) = - \int_0^{\omega_c} \mathrm{d} \omega J(\omega)\frac{(1 - \cos(\omega t))}{\omega^2} \coth(β\omega/2) ,
```
with ``\beta = (k_B T)^{-1}``. For the case where ``\beta \longrightarrow \infty``, the integral reads
```math
    \Gamma(t) = - \int_0^{\omega_c} \mathrm{d} \omega J(\omega)\frac{(1 - \cos(\omega t))}{\omega^2} ,
```
The time-dependent anti-diagonal terms of the reduced density matrix are then expressed as

```math
    \rho_{12}(t) = \rho_{21}(t)^* =\rho_{12}(0) \exp(\Gamma(t)) 
```

Setting up the initial two-level system as a cat state ``|\psi\rangle_S(0) = \frac{|0\rangle \pm |1\rangle}{\sqrt{2}}``, this leads to ``\rho_{12}(0)=\frac{1}{2}``.

Here we break out and comment the script in `MPSDynamics/examples/puredephasing.jl` and `MPSDynamics/examples/puredephasing_temperature.jl` to show how to simulate this model with an Ohmic SD (hard cut-off) using the T-TEDOPA method as implemented in `MPSDynamics.jl`.

The T-TEDOPA method relies on a truncated chain mapping that transform the initial Hamiltonian into
```math
        \hat{H} = \frac{\omega_0}{2} \hat{\sigma}_z +  c_0 \frac{\hat{\sigma}_z}{2}(\hat{b}_0^\dagger + \hat{b}_0) + \sum_{i=0}^{N-1} t_i (\hat{b}_{i+1}^\dagger \hat{b}_i + \mathrm{h.c.}) + \sum_{i=0}^{N-1} \epsilon_i \hat{b}_i^\dagger \hat{b}_i
```

## The code

First, we load the `MPSdynamics.jl` package to be able to perform the simulation, the `Plots.jl` one to plot the results, the `LaTeXStrings.jl` one to be able to use ``\LaTeX`` in the plots and eventually `QuadGK.jl` to perform the analytical integral calculations.

```julia
using MPSDynamics, Plots, LaTeXStrings, QuadGK
```
We then define variables for the physical parameters of the simulation.
Among these, two are convergence parameters:

*  `d` is the number of states we retain for the truncated harmonic oscillators representation of environmental modes
* `N` is the number of chain (environmental) modes we keep. This parameters determines the maximum simulation time of the simulation: indeed excitations that arrive at the end of the chain are reflected towards the system and can lead to unphysical results

```julia
#----------------------------
# Physical parameters
#----------------------------

ΔE = 0.008 # Energy of the electronic states

d = 5 # number of Fock states of the chain modes

N = 30 # length of the chain

α = 0.01 # coupling strength

s = 1 # ohmicity

ωc = 0.035 # Cut-off of the spectral density J(ω)

```
The chain coefficient have still to be calculated. This part is the only difference in parameters for the zero-temperature case and the thermalized bath. For ``T = 0 ~ \text{K}``, the chain coeffficients can be calculated with the function [`MPSDynamics.chaincoeffs_ohmic`](@ref) :
```julia
cpars = chaincoeffs_ohmic(N, α, s; ωc=ωc) # chain parameters, i.e. on-site energies ϵ_i, hopping energies t_i, and system-chain coupling c_0
``` 
Concerning the case ``T \neq 0 ~ \text{K}``, after having set up a temperature value, the chain coefficients can be calculated with the function [`MPSDynamics.chaincoeffs_finiteT`](@ref)

```julia
β = 100

cpars = chaincoeffs_finiteT(N, β; α=α, s=s, J=nothing, ωc=ωc, mc=4, mp=0, AB=nothing, iq=1, idelta=2, procedure=:Lanczos, Mmax=5000, save=false)  # chain parameters, i.e. on-site energies ϵ_i, hopping energies t_i, and system-chain coupling c_0
```

We set the simulation parameters and choose a time evolution method.
As always for simulations of dynamics, the time step must be chosen wisely. The error of the TDVP methods is ``\mathcal{O}(dt^3)``.
In this example we present two one-site implementation of TDVP that both preserves the unitarity of the evolution:

* the regular one-site method with the keyword `:TDVP1` where all the virtual bonds of the MPS have the same bond dimension ``D``
* the adaptive method with the keyword `:DTDVP` where the bond dimension is locally increased at each time step if the TDVP projection error crosses a threshold value

Logically the constant bond dimension of the MPS for TDVP1 and the threshold of the projection error for DTDVP are their respective convergence parameter.
```julia
#-----------------------
# Simulation parameters
#-----------------------

dt = 1.0 # time step

tfinal = 300.0 # simulation time

method = :TDVP1 # time-evolution method

D = 2 # MPS bond dimension
```
Using `MPSDynamics.jl` built-in methods we define the Pure Dephasing model MPO and the MPS representing the initial state.
This initial state is a product state between the system and the chain. It is constructed using a list of the 'local state' of each site of the MPS, and the dimensions of the physical legs of the MPS are set to be the same as the ones of the MPO.
```julia
#---------------------------
# MPO and initial state MPS
#---------------------------

H = puredephasingmpo(ΔE, d, N, cpars)

# Initial electronic system in a superposition of 1 and 2
ψ = zeros(2)
ψ[1] = 1/sqrt(2)
ψ[2] = 1/sqrt(2)

A = productstatemps(physdims(H), state=[ψ, fill(unitcol(1,d), N)...]) # MPS representation of |ψ>|Vacuum>
```
We then chose the observables that will be stored in the data and the [`MPSDynamics.runsim`](@ref) arguments. This example relies on the storage of the reduced density matrix (called `reduceddensity` in runsim).
```julia

#---------------------------
# Definition of observables
#---------------------------

ob1 = OneSiteObservable("sz", sz, 1)


#-------------
# Simulation
#------------

A, dat = runsim(dt, tfinal, A, H;
                name = "pure dephasing model with temperature",
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
```
After having performed the dynamics, the analytical decoherence function is numerically calculated with the help of the quadgk function (from the `QuadGK.jl` package). For the case ``T = 0 ~ \text{K}`` it reads
```julia
#----------
# Analytical results at specified temperature 
# (see The Theory of Open Quantum System, H.-P. Breuer & F. Petruccione 2002, Chapter 4)
#----------

Johmic(ω,s) = (2*α*ω^s)/(ωc^(s-1))

time_analytical = LinRange(0.0,tfinal,Int(tfinal))

Γohmic(t) = - quadgk(x -> Johmic(x,s)*(1 - cos(x*t))/x^2, 0, ωc)[1]

Decoherence_ohmic(t) = 0.5*exp(Γohmic(t))
```
whereas the analytical integral is changed for the case ``T \neq 0 ~ \text{K}``
```julia
Γohmic(t) = - quadgk(x -> Johmic(x,s)*(1 - cos(x*t))*coth(β*x/2)/x^2, 0, ωc)[1]
```

Eventually, the stored reduced density matrix is compared against the analytical formula
```julia
#-------------
# Plots
#------------

ρ12 = abs.(dat["data/Reduced ρ"][1,2,:])

plot(time_analytical, t->Decoherence_ohmic(t), label="Analytics", title=L"Pure Dephasing, Ohmic $s=%$s$, $\beta = %$β ~\mathrm{K}$", linecolor=:black, xlabel="Time (arb. units)", ylabel=L"Coherence $|\rho_{12}(t)|$", linewidth=4, titlefontsize=16, legend=:best, legendfontsize=16, xguidefontsize=16, yguidefontsize=16, tickfontsize=10)

plot!(dat["data/times"], ρ12, lw=4, ls=:dash, label="Numerics")
```

## Bibliography

[^breuer]:
      > Breuer, H.;Petruccione, F. The Theory of Open Quantum Systems; Oxford University Press, 2002.


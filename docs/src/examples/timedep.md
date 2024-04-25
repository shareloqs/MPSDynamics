# Time-dependent Hamiltonian

## Context

To simulate a drive or a laser pulse, a time-dependent Hamiltonian is often needed. Here we explain thanks to the script in `MPSDynamics/examples/sbm_Htimedependent.jl` how this type of Hamiltonian can be integrated within the method and applied on a wavefunction. We will take a Spin-Boson Model (SBM) Hamiltonian with an Ohmic spectral density. However, the time-dependency is not model specific and can be adapated for other models. For more information about the SBM, see the dedicated example. For this example, we simulate a drive of the form

```math
        \hat{H}_\text{drive}(t) = \epsilon \hat{\sigma}_x \sin(\omega_\text{drive} t)
```
with ``\epsilon = \frac{2 \pi}{T_\text{Rabi}}`` the amplitude of the drive, ``T_\text{Rabi}`` the Rabi period, and ``\omega_\text{drive}`` the frequency of the drive. The drive is set up to be on resonance with the two-level system.

## The code
First we load the `MPSdynamics.jl` package to be able to perform the simulation, the `Plots.jl` one to plot the results, and the `LaTeXStrings.jl` one to be able to use ``\LaTeX`` in the plots. The function [`MPSDynamics.disp`](@ref) is also imported.

```julia
using MPSDynamics, Plots, LaTeXStrings

import MPSDynamics: disp
```

We then define variables for the physical parameters of the simulation.
Among these, three are convergence parameters:

*  `d` is the number of states we retain for the truncated harmonic oscillators representation of environmental modes
* `N` is the number of chain (environmental) modes we keep. This parameters determines the maximum simulation time of the simulation: indeed excitations that arrive at the end of the chain are reflected towards the system and can lead to unphysical results

The variable `Ndrive` represents the site of the MPO where the operator of the time-dependent part acts. For this example, the two-level system is at the first site of the MPS. 

```julia
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
```
We set the simulation parameters and choose a time evolution method.
As always for simulations of dynamics, the time step must be chosen wisely. The error of the TDVP methods is ``\mathcal{O}(dt^3)``.
In this example we present only one-site implementation of TDVP that preserves the unitarity of the evolution:

* the regular one-site method with the keyword `:TDVP1` where all the virtual bonds of the MPS have the same bond dimension ``D``

Logically the constant bond dimension of the MPS for TDVP1 is the respective convergence parameter.

```julia
#-----------------------
# Simulation parameters
#-----------------------

dt = 0.5 # time step

tfinal = 100.0 # simulation time

method = :TDVP1 # time-evolution method

D = [6] # MPS bond dimension
```
Using `MPSDynamics.jl` built-in methods we define the SBM MPO and the MPS representing the initial state.
This initial state is a product state between the system and the chain. It is constructed using a list of the 'local state' of each site of the MPS, and the dimensions of the physical legs of the MPS are set to be the same as the ones of the MPO.


In this part, the time-dependent terms of the MPO are stored in a list. This enables to add these terms to the MPO elements at each timestep and avoid the calculation of an entire new MPO. This way is not cumbersome since it adds a matrix sum at each time step.

```julia
#---------------------------
# MPO and initial state MPS
#---------------------------

timelist = collect(0:dt:tfinal)
numsteps = length(timelist)-1

Ht = [ϵ*sx*sin(ωdrive*tstep) for tstep in timelist] # Time-dependent Hamiltonian term

H = spinbosonmpo(ω0, Δ, d, N, cpars) # MPO representation of the Hamiltonian

ψ = unitcol(2,2) # Initial down-z system state 

A = productstatemps(physdims(H), state=[ψ, fill(unitcol(1,d), N)...]) # MPS representation of |ψ>|Vacuum>
```
We then choose the observables that will be stored in the data and the [`MPSDynamics.runsim`](@ref) arguments.
```julia
#---------------------------
# Definition of observables
#---------------------------

ob1 = OneSiteObservable("sz", sz, 1)

ob2 = OneSiteObservable("chain mode occupation", numb(d), (2,N+1))

ob3 = TwoSiteObservable("SXdisp", sx, disp(d), [1], collect(2:N+1))

```
[`MPSDynamics.runsim`](@ref) is called to perform the dynamics. The argument `timedep` is set up to true and `Ndrive` and `Htime` are provided in the kwargs. 
```julia
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
                timedep = true, # the Hamiltonian is time dependent
                Ndrive = Ndrive, # the first site of the MPS/MPO (i.e. the system) is concerned
                Htime = Ht, # list of time-dependent terms
                verbose = false,
                save = true,
                plot = true,
                );
```
Eventually, the stored observables can be represented
```julia

#----------
# Plots
#----------

plot(dat["data/times"], dat["data/sz"], label="Dmax = $(D...)", xlabel=L"t",ylabel=L"\sigma_z", title="")
```


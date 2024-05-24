# User Guide

Here we explain the different steps to perform a simulation.

Examples with detailed explanations can be found in [Examples](@ref).

## Initial State

The initial many-body state of the {System + Environment} can be described easily as [Matrix Product States](@ref) (MPS) or [Tree Tensor Networks](@ref) (TTN) state (e.g. when a system is coupled to several environments).

### Matrix Product States

A MPS can be initialized with several methods.

The [`productstatemps`](@ref) method enables  to instantiate arbitrary MPS of fixed uniform bond dimension with non-uniform physical dimensions.
The individual states of the MPS sites can be provided by setting state to a list of column vectors. 
Setting `state=:Vacuum` will produce an MPS in the vacuum state. 
Setting `state=:FullOccupy` will produce an MPS in which each site is fully occupied.
The gauge of the MPS can also be set using a keyword argument.

```julia
julia> ψ = unitcol(1,2) # system initial state

julia> d = 6; N = 30; α = 0.1; Δ = 0.0; ω0 = 0.2; s = 1 

julia> cpars = chaincoeffs_ohmic(N, α, s) # chain coefficient for an Ohmic spectral density

julia> H = spinbosonmpo(ω0, Δ, d, N, cpars)

julia> A = productstatemps(physdims(H), state=[ψ, fill(unitcol(1,d), N)...]) # MPS representation of |ψ>|Vacuum>
```

Alternatively, a chain with a specified number of excitation localised on one site, or delocalized accross several sites can be generated with [`MPSDynamics.chainmps`](@ref).

Random MPS can also be generated with the [`randmps`](@ref) method.

For the case of fermionic states (which need to be anti-symmetrized), the [`MPSDynamics.electronkmps`](@ref) method generate an MPS for an electron with momentum `k`, and the [`MPSDynamics.electron2kmps`](@ref) generate an MPS with 2 electrons in k-states `k1` and `k2`.

### Tree Tensor Networks

Write a quick explanation of the how trees are structured: parents, child, nodes; and that most methods to initialize a MPS are overloaded for `TreeNetwork`. 

## Hamiltonian

In order to perform time evolution and have access to the dynamics of the many-body state a Hamiltonian needs to be specified in the form of a Matrix Product Operator (MPO) of as a tree tensor network.
Either way, this can be done by using a [Build-in Hamiltonian](@ref), [Convert a MPO from ITensor](@ref), or creating a [Tailored MPO](@ref).

In the context of Open Quantum Systems, custom chain coefficients for the environment can be generated for finite temperature simulations, and/or user provided spectral densities (SDs).

### Build-in Hamiltonian

MPSDynamics provides several topical Hamiltonians directly in the form of MPO or Tree Tensor Networks such as the Ising model [`MPSDynamics.isingmpo`](@ref), the XYZ Hamiltonian [`MPSDynamics.xyzmpo`](@ref), the Spin Boson Model [`MPSDynamics.spinbosonmpo`](@ref), a spin coupled to two bosonic baths `MPSDynamics.twobathspinmpo`, nearest neighbour interactions Hamiltonian `MPSDynamics.nearestneighbourmpo`, the independent boson model `MPSDynamics.ibmmpo`, (non-)uniform tight-binding chain Hamiltonian [`MPSDynamics.tightbindingmpo`](@ref).

### Convert a MPO from ITensor

The method [`MPSDynamics.MPOtoVector`](@ref) converts an ITensors chain MPO into a form compatible with MPSDynamics.

### Tailored MPO

One can also construct MPO tailored to the problem one is interested in.
MPOs are fundamentally a lists of rank-4 tensors such that the right bond dimension of the nth tensor must be equal to the left bond dimension of the n+1th tensor; and the dimension of the physical bonds of the nth tensor must be equal to the corresponding physical bond on the MPS.

### Finite Temperature and Custom SD

In the T-TEDOPA framework (see [Theoretical Background](@ref) for more details) finite temperature simulations are done with an effective pure state description of the system and the environment where the coupling coefficients (or the SD) is temperature-dependent.

The corresponding chain coefficients for an Ohmic or a user provided spectral density (that can thus in pratice be either at zero or finite temperature) are computed with the `[`chaincoeffs_finiteT`](@ref)`.
This method is based on the `ORTHOPOL` routines[^Gautschi]

## Observables

System and environment observables can be computed, as well as system-and-environment 'non-local' observables.

Observables that will be passed to `MPSDynamics.runsim`(@ref) to have their expectation value computated at each time step are defined with the [`OneSiteObservable`](@ref) and `[`TwoSiteObservable`](@ref)`.

One-site and two-site obsevables work similarly, they need to be given a name, an (pair of) operator(s) and the (list of) site(s) on which they are evaluated.

For instance one can calculated the average number of excitation in each of the ``N`` environmental modes

```julia
    ob = OneSiteObservable("chain mode occupation", numb(d), (2,N+1))
```

It is also possible to measure composite system/environment observables, for example

```julia
    ob = TwoSiteObservable("SXdisp", sx, disp(d), [1], collect(2:N+1))
```
which measure the correlation between the spin in the x-direction and the displacement of the bath modes.

Purely environmental 'non-local' observables such as 

```julia
    ob = TwoSiteObservable("bath coherence", crea(d), anih(d), collect(2:N+1), collect(2:N+1))
```
that computes all the chain mode coherences ``\langle\hat{a}_n^\dagger\hat{a}_{m}\rangle`` (and the chain mode occupation when ``n = m``).

We note that if one knows the coherences and populations of all the chain modes, it is then possible to reconstruct the populations of the normal mode environment.

## Time-Evolution

Simulation are performed when calling the [`MPSDynamics.runsim`](@ref) function where the time-evolution method should be specified.

The time-evolution methods currently implemented belong to the familly of Time-Dependent Variational Principle (TDVP).
The central point of this method, in the modern tensor networks formulation, is that instead of solving the Schrödinger equation and then truncating the MPS representation of the quantum state, one can solve the equations of motion projected into a space of restricted bond dimension. 
The major advantage of this method is that it naturally preserves the unitarity of the time evolution and conserves the energy (except in its two-site implementation).
Three variants of TDVP are implemented in the `MPSDynamics.jl` package:
* one-site TDVP (TDVP1): fixed bond dimension, preserves unitarity, conserves energy, scales as ``\mathcal{O}(D^2 d^2 w^2 + D^3 dw + D^3 d^2 )`` where ``D`` is the MPS (or TTN) bond dimension, ``d`` the local dimension, ``w`` the MPO bond dimension.
* two-site TDVP (TDVP2): adaptive bond dimension, breaks unitarity, scales as ``\mathcal{O}(D^2d^3 w^2 + D^3 d^2 w + D^3 d^3 )``.
* adaptive one-site TDVP (DTDVP): adaptive bond dimension, preserves unitarity, conserves energy.

Local one-site and two-site observables, as well as non-local two-site observables, can be efficiently computed for each time-step of the time-evolution with the three methods.

### One-site TDVP

This implementation of the TDVP can be used for MPS and TTN.
It performs one-site updates of the tensor network tensors while keeping the bond dimension fixed.
Because no reduction of the dimensionality is performed (for instnce by throwing away some simgular values), it naturally preserves unitarity and hence conserves energy.
The otherside of the coin is that simulation might be slower because the bond dimension is large at *all time*.

The convergence parameter of the simulation is thus the bond dimension `D`.

This method can be used with in [`MPSDynamics.runsim`](@ref) with the key word argument `method=:TDVP1`.

### Two-Site TDVP

This version of the TDVP evolves two-site tensors than can be split back into two individual site tensors by applying an SVD.
The bond dimension between two neighbouring sites can then be truncated to any value ``≤ r`` by throwing away any singular values that fall below some threshold. 
In this way, the MPS bond dimension will grow dynamically throughout the course of the evolution to capture entanglement, as and when it emerges.
The truncation entails a loss of unitarity. 
Indeed, for a nearest-neighbour Hamiltonian, applying the two-site projector will not entail a projection error, leading to a scheme that is almost identical to TEBD,wherein the error arises solely from the truncation.
Furthermore, 2TDVP scales poorly with the local physical dimension, and is known to have issues with long-range interactions.

The convergence parameter is the threshold of the SVD.

This method can be used with in [`MPSDynamics.runsim`](@ref) with the key word argument `method=:TDVP2`.

### Adaptive One-Site TDVP

Ahead of any time evolution, one computes new bond dimensions of the MPS if the relative rate of change of the TDVP projection error with respect to the bond dimension is larger than a chosen precision `p`.
Then, TDVP1 is performed, using projectors with sub-spaces expanded accordingly, to produce an MPS evolved by one time step with the new,increased, bond dimensions.
This version is faster than TDVP1 due to the acceleration gained from having more optimised bond dimensions; the bond update step is cheap and so its cost should not normally outweigh this advantage.
However, two things need to be noticed:
* DTDVP can get stuck into local minima from the initial time. If the bond dimension does not change during the time-evolution, consider embeding your MPS in a larger manifold with [`mpsembed!`](@ref) before time-evolving it.
* For models with analytical solutions that are MPS of a given bond dimension (such as the independent boson model), DTDVP can overshoot the analytical bond dimension because the relative rate of change of the projection error becomes dominated by random numerical fluctuations.

The convergence parameter is a threshold value `p` for the rate of change of the projection error with respect to the bond dimension.

This method can be used with in [`MPSDynamics.runsim`](@ref) with the key word argument `method=:DTDVP`.

## Data Storage

The data (i.e. observables time-series) is stored in the JLD format which is based on HDF5.

The HDF5 format is natively supported across many platforms and languages (e.g. `Python`, or `Mathematica`).

For the data to be saved to a file after a run, the keyword argument `save=true` needs to be used in `runsim`.

The directory where the data should be saved can be chosen by setting a path with the `savedir` keyword argument.

A `plot` keyword argument can also be used to choose whether plots for 1D observables will be automatically generated and saved along with the data. 

Loading the data in Julia using the [`JLD.jl`](https://github.com/JuliaIO/JLD.jl) package will recover the full type information of the Julia variables that were stored.

```julia
julia> using JLD

julia> dat = load("filename.jld")
Dict{String, Any} with 2 entries:
  "parameters" => Dict{String, Any}("tmax"=>0.2, "method"=>:DTDVP, "dt"=>0.0005…
  "data"       => Dict{String, Any}("bonddims"=>[1 1 … 1 1; 1 2 … 2 2; … ; 1 1 …
```

If the data is loaded in an variable `dat` the structure is the folowing:
- the `parameters` entry contains a dictionary where all the simulation parameters are stored
- the `data` entry contains a dictionary where are stored simulation time, the observables and (whenever relevent) the bond dimension of the state at each time steps. 

```julia
julia> dat["parameters"]
Dict{String, Any} with 11 entries:
  "tmax"       => 0.2
  "method"     => :DTDVP
  "dt"         => 0.0005
  "name"       => "my model"
  "Δ"          => 0.0
  "β"          => 0.0186854
  "N"          => 300.0
  "d"          => 15.0
  "unid"       => "Ovzm6"
  "ω0"         => 0.0
  "convparams" => 0.0005

julia> dat["data"]
Dict{String, Any} with 6 entries:
  "bonddims" => [1 1 … 1 1; 1 2 … 2 2; … ; 1 1 … 7 7; 1 1 … 1 1]
  "sx"       => [1.0, 0.912818, 0.741759, 0.605797, 0.528792, 0.492497, 0.47976…
  "sz"       => [0.0, 0.0871825, 0.25824, 0.394201, 0.471207, 0.507503, 0.52023…
  "nchain"   => [0.0 0.23466 … 1.84319 1.76098; 0.0 0.00231507 … 0.83105 0.9033…
  "sy"       => [0.0, -0.0133489, -0.0588887, -0.0858181, -0.0759996, -0.048539…
  "times"    => [0.0, 0.0005, 0.001, 0.0015, 0.002, 0.0025, 0.003, 0.0035, 0.00…

```

[^Gautschi]:
	> Gautschi, W. Algorithm 726: ORTHPOL–a package of routines for generating orthogonal polynomials and Gauss-type quadrature rules. ACM Trans. Math. Softw. 20, 21–62 (1994).


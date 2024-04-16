# User Guide

Here we explain the different steps to perform a simulation.

Examples with detailed explanations can be found in .

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

### Build-in Hamiltonian

MPSDynamics provides several topical Hamiltonians directly in the form of MPO or Tree Tensor Networks such as the Ising model [`MPSDynamics.isingmpo`](@ref), the XYZ Hamiltonian [`MPSDynamics.xyzmpo`](@ref), the Spin Boson Model [`MPSDynamics.spinbosonmpo`](@ref), a spin coupled to two bosonic baths `MPSDynamics.twobathspinmpo`, nearest neighbour interactions Hamiltonian `MPSDynamics.nearestneighbourmpo`, the independent boson model `MPSDynamics.ibmmpo`, (non-)uniform tight-binding chain Hamiltonian [`MPSDynamics.tightbindingmpo`](@ref).

### Convert a MPO from ITensor

The method [`MPSDynamics.MPOtoVector`](@ref) converts an ITensors chain MPO into a form compatible with MPSDynamics.

### Tailored MPO

One can also construct MPO tailored to the problem one is interested in.
MPOs are fundamentally a lists of rank-4 tensors such that the right bond dimension of the nth tensor must be equal to the left bond dimension of the n+1th tensor; and the dimension of the physical bonds of the nth tensor must be equal to the corresponding physical bond on the MPS.

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

Explain that we do TDVP.

A simulation is runned by the `MPSDynamics.runsim`(@ref) function.

Local one-site and two-site observables, as well as non-local two-site observables, can be efficiently computed for each time-step of the time-evolution.


### One-site TDVP

Fixed bond dimension, complexity XXX, preserves unitarity.

The convergence parameter is the bond dimension.

Is set in `MPSDynamics.runsim`(@ref) usind the key word argument `method=:TDVP1`

### Two-Site TDVP

Varying bond dimension, complexity higher than 1TDVP, but breaks unitarity because of a SVD.

The convergence parameter is the threshold of the SVD.

Is set in `MPSDynamics.runsim`(@ref) usind the key word argument `method=:TDVP2`

### Adaptive One-Site TDVP

Keeps the good scaling of the 1TDVP, preserves unitarity but is able to increase the bond dimendion.

The convergence parameter is a threshold value for the rate of change of the projection error with respect to the bond dimension.

Is set in `MPSDynamics.runsim`(@ref) usind the key word argument `method=:DTDVP`

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

# Users Guide

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
julia> ψ = unitcol(1,2); d = 6; N = 30; α = 0.1; Δ = 0.0; ω0 = 0.2; s = 1

julia> cpars = chaincoeffs_ohmic(N, α, s)

julia> H = spinbosonmpo(ω0, Δ, d, N, cpars)

julia> A = productstatemps(physdims(H), state=[ψ, fill(unitcol(1,d), N)...]) # MPS representation of |ψ>|Vacuum>
```

Alternatively, a chain with a specified number of excitation localiswed on one site, or delocalized accross several sites can be generated with [`MPSDynamics.chainmps`](@ref).

Random MPS can also be generated with the [`randmps`](@ref) method.

For the case of fermionic states (which need to be anti-symmetrized), the [`MPSDynamics.electronkmps`](@ref) method generate an MPS for an electron with momentum `k`, and the [`MPSDynamics.electron2kmps`](@ref) generate an MPS with 2 electrons in k-states `k1` and `k2`.

### Tree Tensor Networks

Write a quick explanation of the how trees are structured: parents, child, nodes; and that most methods to initialize a MPS are overloaded for `TreeNetwork`. 

## Hamiltonian

In order to perform time evolution and have access to the dynamics of the many-body state a Hamiltonian needs to be specified in the form of a Matrix Product Operator (MPO) of as a tree tensor network.
Either way, this can be done by using a [Build-in Hamiltonian](@ref), [Convert a MPO from ITensor](@ref), or creating a [Tailored MPO](@ref).

### Build-in Hamiltonian

MPSDynamics provides several topical Hamiltonians directly in the form of MPO or Tree Tensor Networks such as the Ising model [`MPSDynamics.isingmpo`](@ref), the XYZ Hamiltonian [`MPSDynamics.xyzmpo`](@ref), the Spin Boson Model [`MPSDynamics.spinbosonmpo`](@ref), (non-)uniform tight-binding chain Hamiltonian [`MPSDynamics.hbathchain`](@ref).

### Convert a MPO from ITensor

The method [`MPSDynamics.MPOtoVector`](@ref) converts an ITensors chain MPO into a form compatible with MPSDynamics.

### Tailored MPO

Explain that a MPO is just a list of tensors with matching bond dimensions.

## Observables

## Time-Evolution

## Data Storage

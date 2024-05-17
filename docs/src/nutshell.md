# In a nutshell: MPSDynamics.jl for open quantum systems


`MPSDynamics.jl` was originally developed to perform tensor network simulations of open quantum systems by solving the Schrödinger equation for the closed {System + Environment}. 
The key idea to construct the simulations is to reformulate the open quantum system Hamiltonian as a one-dimensional many-body problem with nearest-neighbor interactions. 
The dynamics are then efficiently simulated with [Tensor Networks](@ref) methods. 

![Sketch of the different ingredients behind MPSDynamics](examples/mappings-1.png)

The starting point is a system linearly coupled to a bosonic environment, where the coupling between system and environment is specified by the spectral density function $J(\omega)$: the prototypical example being [The Spin-Boson Model](@ref).
At finite temperature, the initial state of the environment will be a thermal state (mixed), requiring density matrix formalism and thus making the problem more complex.
To circumvent this issue, we consider instead an _extended environment_, with coupling to the system characterized by a thermalized spectral density function $J(\omega,\beta)$. 
The two environments induce the same reduced dynamics on the system, enabling us to deal with pure states only[^1].
The vacuum state (pure) of the extended environment corresponds to the thermal state of the original environment. 

To exploit the computational efficiency of the matrix product states (MPS) description of pure states, we apply the so called [Chain-Mapping of bosonic environments](@ref) —a unitary transformation dependent on $J(\omega,\beta)$— to the Spin-Boson Hamiltonian, mapping it on a chain Hamiltonian with nearest-neighbor interactions.
This enables efficient representation of the full {System + Environment} state and its time evolution using the [Time-Dependent Variational Principle](@ref).

To go further, you will find in this documentation:
* a [User Guide](@ref) explaining how to use the package
* several examples illustrating different features of the package: adaptive time-evolution methods, [Time-dependent Hamiltonian](@ref), or fermionic environments with [The Anderson Impurity Model](@ref)
* a [Theoretical Background](@ref) covering the chain mapping procedure, flying over the basics of tensor networks, and summarizing the TDVP methods

If this package was useful in your work, do not forget [Citation](@ref).
And if you would like to get involved in its development, you can find out [How to Contribute](@ref).

[^1]: Moreover, the dynamics of the original environment can still be recovered, see [Inspecting the bath by undoing the chain mapping](@ref)

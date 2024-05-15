# Inspecting the bath by undoing the chain mapping

Here we give some context on the example script provided in `MPSDynamics/example/bath-observables.jl`. This example demonstrates the setup and execution of a simulation for a two-level system coupled to an Ohmic bath at finite temperature, where we exploit the access to the chain observables to:
- undo the chain mapping[^chin_exact_2010], thus obtaining their representation in the extended bath of T-TEDOPA [^tamascelli_efficient_2019], characterized by $J(\omega, \beta)$
- inverting the thermofield transformation[^devega_thermo_2015], thus obtaining the representation of the physical frequencies in the original environment, characterized by $J(\omega)$

T-TEDOPA allows to substitute a thermally occupied bath by an extended one, in the pure state of the vacuum, extending the bath of frequencies to negative values: the _creation_ from the system of a mode of _negative frequency_ in the extended bath of frequencies corresponds to the _absorption_ of energy for the system from the thermally occupied modes of the environment. However, by exploiting the thermofield transformation, it is possible to recover the picture of the occupations of the steady states in the physical bath of frequencies[^riva_thermal_2023]. 

## The code

We start by defining the parameters of the simulation:

```julia
d = 10      # number of Fock states of the chain modes
N = 60      # length of the chain
α = 0.01    # coupling strength
ω0 = 0.2  # TLS gap
s = 1       # ohmicity
ωc = 1.  # Cut-off of the spectral density J(ω)
β = 20    # Thermalized environment
```
We set the specifics of the simulation:
```julia
method = :TDVP1         # time-evolution method
conv = 3                # bond dimension for the TDVP1
dt = 0.5                # time step
tfinal = 60.0           # simulation time
Tsteps = Int(tfinal / dt)
```
And then compute the chain coefficients, i.e. on-site energies $\epsilon_i$, hopping energies $t_i$, and system-chain coupling $c_0$, that define the chain representation of the pure-dephasing model:
```julia
cpars = chaincoeffs_finiteT(N, β; α=α, s=s, J=nothing, ωc=ωc, mc=4, mp=0, AB=nothing, iq=1, idelta=2, procedure=:Lanczos, Mmax=5000, save=false)
```
with this, we can compute the corresponding MPO, and the initial state in MPS form ($1/\sqrt{2}(|0\rangle + |1\rangle)$):
```julia
H = puredephasingmpo(ω0, d, N, cpars)

ψ = zeros(2)
ψ[1] = 1/sqrt(2)
ψ[2] = 1/sqrt(2)

A = productstatemps(physdims(H), state=[ψ, fill(unitcol(1,d), N)...])
```
We can now define the observables we are interested in computing. Importantly, we not only compute observables related to the system, but also to the chain modes, so that we can inspect the environment.
```julia
ob1 = OneSiteObservable("sz", sz, 1)
ob2 = OneSiteObservable("sx", sx, 1)
ob3 = OneSiteObservable("chain_mode_occupation", numb(d), (2,N+1))
ob4 = OneSiteObservable("c", crea(d), collect(2:N+1))
ob5 = OneSiteObservable("cdag", crea(d), collect(2:N+1))
ob6 = TwoSiteObservable("cdagc", crea(d), anih(d), collect(2:N+1), collect(2:N+1))
ob7 = TwoSiteObservable("cdagcdag", crea(d), crea(d), collect(2:N+1), collect(2:N+1))
ob8 = TwoSiteObservable("cc", anih(d), anih(d), collect(2:N+1), collect(2:N+1))
```
with this we run the simulation:
```julia

```

Furthermore, in the pure dephasing case, it is also possible to obtain the analytical prediction of the time evolution of the occupations of the bath's modes, so that we can compare our numerical results with the analytical ones, exploiting the Heisenberg time evolution relation:
$$
\frac{d \langle \hat b_\omega \rangle}{dt} = -i \langle[ \hat b_\omega, \hat H] \rangle = - i \omega \langle\hat b_\omega \rangle - i \frac{\langle \hat \sigma_z \rangle}{2} \sqrt{J(\omega, \beta)},     \\
\frac{d \langle \hat n_\omega \rangle}{dt} = -i \langle[\hat b_\omega^\dagger \hat b_\omega, \hat H] \rangle= 2 \frac{|J(\omega,\beta)|}{\omega} \sin(\omega t).
$$
To this end, it is convenient to choose one of the eigenstates of $\hat \sigma_z$ as the initial state, so that $\langle \hat \sigma_z \rangle = \pm 1$. By solving these differential equations, one obtains the time evolved theoretical behavior of the bath.

During a numerical simulation we work with chain of finite length $N$. This truncation on chain modes, introduces a sampling on the modes in the original star-like environment. To recover the frequency modes that are implicitly sampled, one has to diagonalize the tri-diagonal $N\times N$ matrix $H^\text{chain}$, where the diagonal is formed by the $e_n$ coefficients, that is the chain's frequencies, and the upper and lower diagonals by the $N-1$ hopping coefficients $t_n$. The unitary matrix that defines the change of basis from the star-like to the chain-like environment is $U_n$, constituted by the eigenvectors of $H^\text{chain}$. The reverse path could also be taken: one can decide to sample the relevant frequencies in the star-like environment, and apply the inverse transformation to construct the tri-diagonal matrix defining the chain. Once the problem is mapped on the chain, the MPO representation of the new Hamiltonian follows straightforwardly.




___________________
# References

[^chin_exact_2010]:
    > Chin, A. W.; Rivas, Á.; Huelga, S. F.; Plenio, M. B. Exact Mapping between System-Reservoir Quantum Models and Semi-Infinite Discrete Chains Using Orthogonal Polynomials. Journal of Mathematical Physics 2010, 51 (9), 092109. https://doi.org/10.1063/1.3490188.

[^tamascelli_efficient_2019]:
    > Tamascelli, D.; Smirne, A.; Lim, J.; Huelga, S. F.; Plenio, M. B. Efficient Simulation of Finite-Temperature Open Quantum Systems. Phys. Rev. Lett. 2019, 123 (9), 090402. https://doi.org/10.1103/PhysRevLett.123.090402.

[^devega_thermo_2015]:
    > de Vega, I.; Banuls, M-.C. Thermofield-based chain-mapping approach for open quantum systems. Phys. Rev. A 2015, 92 (5), 052116. https://doi.org/10.1103/PhysRevA.92.052116.

[^riva_thermal_2023]:
    > Riva, A.; Tamascelli, D.; Dunnett, A. J.; Chin, A. W. Thermal cycle and polaron formation in structured bosonic environments. Phys. Rev. B 2023, 108, 195138, https://doi.org/10.1103/PhysRevB.108.195138.

# Theoretical Background

## Chain-Mapping of bosonic environments

We consider, in the Schrödinger picture, a general Hamiltonian where a non-specified system interacts linearly with a bosonic environments

```math
\begin{aligned}
    \hat{H} =& \hat{H}_S + \int_0^{+\infty} \hbar\omega\hat{a}^\dagger_\omega\hat{a}_\omega\ \mathrm{d}\omega + \hat{A}_S\int_0^{+\infty}\sqrt{J(\omega)}\left(\hat{a}_\omega + \hat{a}^\dagger_\omega\right)\mathrm{d}\omega
\end{aligned}
```


where $\hat{a}_\omega$ ($\hat{a}^\dagger_\omega$) is a bosonic annihilation (creation) operator for a normal mode of the environment of energy $\hbar\omega$, $\hat{A}_S$ is a system operator, and $J(\omega) = \sum_k |g_k|^2\delta(\omega - \omega_k)$ is the bath spectral density (SD), defined with the microscopic system-environment coupling strength $g_k$.
The SD quantifies the coupling strengths of the different normal modes of the environment with the system.
Any SD that is not flat corresponds to a non-Markovian environment.

### Zero Temperature

Let us consider the Hamiltonian presented in Eq.(1).
We can introduce a unitary transformation of the continuous normal modes $\hat{a}_\omega$ to an infinite discrete set of interacting modes $\hat{b}_n$[^chin_exact_2010].

```math
    \hat{a}_\omega = \sum_{n=0}^{+\infty} U_n(\omega)\hat{b}_n = \sum_{n=0}^{+\infty} \sqrt{J(\omega)}P_n(\omega)\hat{b}_n\ ,
```

where $P_n(\omega)$ are orthonormal polynomials such that

```math
    \int_{0}^{+\infty}P_n(\omega)P_m(\omega)J(\omega)\mathrm{d}\omega = \delta_{n,m}\ ;
```

and the inverse transformation is

```math
    \hat{b}_n = \int_0^{+\infty} U_n(\omega)\hat{a}_\omega\mathrm{d}\omega\ .
```

Note that the orthonormality of the polynomials ensures the unitarity of the transformation defined in Eq.(2).
The mapping from a continuous set of modes to a (still infinite) discrete set might seem counter-intuitive, however it is a direct consequence of the separability of the underlying Hilbert space.

Under this transformation, the Hamiltonian in Eq.(1) becomes

```math
    \hat{H}= \hat{H}_S + \sum_{n=0}^{+\infty}\varepsilon_n\hat{b}_n^\dagger\hat{b}_n + t_n(\hat{b}_{n+1}^\dagger\hat{b}_n + \mathrm{h.c.}) + \kappa\hat{A}_S(\hat{b}_0 + \hat{b}_0^\dagger)\ .
```

Hence, this mapping transforms the normal bath Hamiltonian into a tight-binding Hamiltonian with on-site energies $\varepsilon_n$ and hopping energies $t_n$.
Another important consequence of this mapping is that now the system only interacts with the first mode $n = 0$ of the chain-mapped environment.
The chain coefficients $\varepsilon_n$, $t_n$, and the coupling $\kappa$ depend solely on the SD.

This makes chain mapping a tool of choice for describing systems coupled to environment with highly structured SD (e.g. experimentally measured or calculated *ab initio*)[^chin_role_2013][^alvertis_nonequilibrium_2019][^dunnett_influence_2021][^caycedosoler_exact_2022].
In this new representation, the Hamiltonian in Eq.(5) has naturally a 1D chain topology.
This makes its representation as a Matrix Product Operator (MPO) and the representation of the joint \{System + Environment\} wave-function as a Matrix Product State (MPS) suited [^orus_practical_2014][^paeckel_timeevolution_2019].

The orthogonal polynomial-based chain mapping and the subsequent representation of the joint wave-function as a MPS (and the operators as MPO) are the building blocks of the Time-dependent Density operator with Orthonormal Polynomials Algorithm (TEDOPA) one of the state-of-the-art numerically exact method to simulate the dynamics of open quantum systems especially in the non-Markovian, non-perturbative regimes both at zero and finite temperatures [^prior_efficient_2010][^woods_simulating_2015][^tamascelli_efficient_2019][^dunnett_simulating_2021][^lacroix_unveiling_2021].

### Finite Temperature with T-TEDOPA

Assuming a unitary evolution for both the system and environment, the system's dynamics can be isolated by tracing out the environmental degrees of freedom. The density operator for the system at time $t$ is described as:

```math
\hat{\rho}_S(t) = \text{T}r_E\left\{\hat{U}(t) \hat{\rho}_S(0) \otimes \hat{\rho}_E(0) \hat{U}^\dagger(t)\right\}.
```

Initially, the system state, $\hat{\rho}_S(0)$, can be pure or mixed, and the environment is in a thermal state defined by the inverse temperature $\beta = (k_B T)^{-1}$. This state is represented by a product of Gaussian states:

```math
\hat{\rho}_E(0) = \bigotimes_\omega \frac{e^{-\beta \omega \hat{b}_\omega^\dagger \hat{b}_\omega}}{Z_\omega(\beta)},
```

The system's evolution is dictated by the environment's two-time correlation function, which in turn is determined by the spectral density function $J$ and the temperature $\beta$:

```math
\hat{S}(t) = \int_0^\infty d\omega J(\omega)\left[e^{-i\omega t}(1 + \hat{n}_\omega(\beta)) + e^{i\omega t} \hat{n}_\omega(\beta)\right],
```

To simulate finite temperature effects using a zero-temperature model with the T-TEDOPA method [^tamascelli_efficient_2019], we extend the spectral density function to cover both positive and negative frequencies, allowing us to use a pure state description for the environment. This is achieved by defining a new spectral density function $J(\omega, \beta)$ that incorporates the Boltzmann factors, supporting the entire real axis:

```math
J(\omega, \beta) = \frac{\text{sign}(\omega)J(\left|\omega\right|)}{2} \Big(1 + \coth\Big(\frac{\beta \omega}{2}\Big)\Big).
```

This modified bath allows us to maintain a pure state description of the environment, represented as a vacuum state, and avoid the computational complexities of density matrices:

```math
\ket{\text{vac}} = \bigotimes_\omega \ket{0}_\omega,
```

The Hamiltonian of the system interacting with this extended bath now includes temperature-dependent interactions:

```math
\hat{H} = \hat{H}_S + \int_{-\infty}^{+\infty} d\omega \omega \hat{b}_\omega^\dagger \hat{b}_\omega + \frac{\hat{\sigma}_x}{2} \otimes \int_{-\infty}^{+\infty} d\omega \sqrt{J(\omega,\beta)}\left(\hat{b}_\omega^\dagger+\hat{b}_\omega\right),
```

This method simplifies the simulation of finite temperature effects by treating them within an effective zero-temperature framework, thereby keeping the computational advantages of using pure states.



## Computation of the chain coefficients

A useful property of the orthonormal polynomials is that they obey a recurrence relation

```math
    P_n(\omega) = (C_{n-1}\omega - A_{n-1})P_{n-1}(\omega) + B_{n-1}P_{n-2}(\omega)\ ,
```

where $A_n$ is related to the first moment of $P_n$, $B_n$ and $C_n$ to the norms of $P_n$ and $P_{n-1}$[^appel_mathematics_2007].
This recurrence relation can be used to construct the polynomials with the conditions that $P_0(\omega) = ||p_0||^{-1} = \left(\int_{\mathbb{R}^{+}} J(\omega)\mathrm{d}\omega \right)^{-\frac{1}{2}}$ and $P_{-1}(\omega) = 0$, with $||\bullet||$ the norm of $\bullet$ with respect to the measure $J(\omega)$, and $P_n(\omega) = p_n(\omega)||p_n||^{-1}$ ; where the polynomials $\{p_n\}_{n\in\mathbb{N}}$ are the so called _monic polynomials_ where the factor $a_n$ in front of $\omega^{n}$ is equal to 1.

The energy of the chain mode $n$ is given by $\varepsilon_n = A_n C_n^{-1}$ and $t_n=C_n^{-1}$ is the coupling between mode $n$ and $n+1$[^chin_exact_2010].

The system couples _only_ to the first mode with the coupling strength $\kappa = ||p_0||$.

Explain that for some weight function/SD they are known analytically and that for others we can use the build-in routines inspired by Gautschi or the PolyChaos.jl package.

## Tensor Networks

A multipartite quantum state $|\psi\rangle$, e.g. a $N$-site system where the sites can each be in a state $|\phi_i\rangle$ belonging to a $d$-dimensional Hilbert space, can be written as follows

```math
    |\psi\rangle = \sum_{\{i_k\}}c_{i_1\ldots i_N}|\phi_{i_1}\rangle\otimes\ldots\otimes|\phi_{i_N}\rangle\ ,
```

where the complex numbers $c_{i_1\ldots i_N}$ are the amplitudes of each state $|\phi_{i_1}\rangle\otimes\ldots\otimes|\phi_{i_N}\rangle$ whose superpositions form in full generality the state $|\psi\rangle$.
Thus the state $|\psi\rangle$ can be completely represented by a rank-$N$ tensor $c$ that is the collection of all possible amplitudes $c_{i_1\ldots i_N}$.
Here by the rank of a tensor, we simply mean the number of indices it has.

### MPS

The tensor $c$ of a quantum state $|\psi\rangle$ corresponding to a one-dimensional system can be decomposed into a product of $N$ smaller rank-3 tensors $T_{k}$ (except for the first and last sites where the tensors will have a rank-2)

```math
    c_{i_1\ldots i_N} = \sum_{\{\alpha\}} T^{\alpha_1}_{i_1}T^{\alpha_1\alpha_2\ }_{i_2}T^{\alpha_2\alpha_3\ }_{i_3}\ldots T^{\alpha_{N-1}}_{i_N} \ .
```

In this form, the local tensor $T_k$ contains the information on the quantum state on site $k$ and its relation (especially the entanglement) with the neighbouring sites.

The decomposition of the tensor of the amplitudes of a quantum state into a product of smaller rank tensors is called a **Matrix Product State** decomposition.

The contracted indices $\alpha_k$ between the tensors are called _virtual indices_ and carry information about the correlations between bi-partitions of the state at bond $k$.
The number of different values a virtual index can take is called the _bond dimension_ and is denoted $D$.
The free indices $i_k$ associated with local quantum states are called _physical indices_.
Thus, they can take $d$ values (with $d$ the dimension of the local Hilbert space).

Any state in the Hilbert space of a one-dimensional many-body system can in principle be represented by a MPS by choosing a sufficiently large value for the bond dimension $D$ \cite{Orus}.
On top of this intellectually satisfying property of MPSs being a dense set of states for a 1d-system, they can also be used as a practical Ansätze for a many-body quantum states by setting a maximal allowed value $\chi$ for the bond dimension $D$.
In doing so, we restrict ourselves to a corner of the total Hilbert space.
The rationale behind this Ansatz is the following: if the initial quantum state of a many-body system has a low bond dimension (typically if the initial state is a product state with $D = 1$), then in a finite time it will only be able to explore a region of the Hilbert space that is not to far away from its starting point.
Thus, the bond dimension will not have the time to diverge exponentially \cite{poulin_quantum_2011}.
However, depending on the physical system at hand, this sub-manifold of the Hilbert space could still be "too large".
There is an additional reason that explains why MPSs are good Ansätze for 1d physical systems.
Most many-body Hamiltonians we (physicists) are interested in are local, meaning that the interactions they describe involve objects that are "neighbours".
For such Hamiltonians, the ground states (outside of potential critical phases) follow the so called _area law_ for the entanglement entropy.\cite{srednicki_entropy_1993, vidal_entanglement_2003, wolf_area_2008}.
This law states that the entanglement entropy $S_{vN}$ of a bi-partition of the system is proportional, not to the volume of the partition as one might expect, but to the hyper-surface of the partition's boundary; hence the name "area law".
For a 3d system this corresponds to an actual surface area $A$, $S_{vN} \sim A$; for a 2d system it corresponds to the length $L$ of the partition's boundary, $S_{vN} \sim L$; and in 1d the boundary reduces to a point, thus the entropy will be independent of the size of the system $S_{vN} \sim \text{constant}$.
The MPSs are states that satisfy this area law.

An application of the [Singular Value Decomposition](https://en.wikipedia.org/wiki/Singular_value_decomposition) is to create efficient approximations of quantum states to perform computations.
The main idea is to reduce the content of the MPS to keep only the parts that contain the physics of interest.
One method to realise this approximation is to do a SVD on each of the tensors of the MPS after each time step of the state time-evolution and to trim the smallest singular values in order to decrease the bond dimension of the MPS down to a chosen maximal value $\chi$.
The corresponding columns and rows of the unitary matrices $U$ and $V^\dagger$ are also removed.
Then, the trimmed matrices $\tilde{U}$, $\tilde{S}$ and $\tilde{V}^\dagger$ are contracted back to give an approximated tensor $T$ with a smaller bond dimension.
Another way to apply the restricted rank approximation is to restrict oneself into working in a manifold of fixed bond dimension $D$ and to use methods that can enforce this constraint.

### MPO

In order to compute expectation values of observables or apply unitary transformations to a quantum state, we need a TN representation of operators.
In the same fashion as a one-dimensional quantum state can be represented as a MPS, operators acting on those states can be represented as **Matrix Product Operators** (MPO).
For an operator $\hat{O}$, its MPO can be defined as follows

```math
    \hat{O} = \sum_{\{i_k\}\{i_k^{'}\} \{w\}} W^{i_1\ i^{'}_1}_{1\ w_0w_1}\ldots  W^{i_N\ i^{'}_N}_{N\ w_{N-1}w_N} |\phi_{i_1^{'}}\ldots \phi_{i_N^{'}}\rangle\langle\phi_{i_1}\ldots \phi_{i_N}| 
```
The contracted indices between the tensors are called _virtual indices_.
The free indices are called _physical indices_ and correspond to the different input and output local quantum states. 
They can take $d$ values (with $d$ the dimension of the local Hilbert space).

### TTN
A natural extension to the MPS is the (loop-free) tree tensor network. 
A TTN is a generalisation of the MPS wherein each site, instead of being connected to only one other site to its right, may be connected to any arbitrary number of _child_ sites.
Provided the tree does not contain any loops, everything that one can do to an MPS/MPO can be extended straight-forwardly to TTN states and TTN operators. 
The generalisation to trees introduces no new conceptual complexity (only implementational complexity).
The sites of a TTN are usually referred to as _nodes_. 
For our purposes, every node of a TTN state and operator has one _parent_ leg, and any number (including zero) of child legs. 
The first node is known as the head-node and has a dummy parent leg with dimension 1.

## Time-Dependent Variational Principal

The original idea behind TDVP goes back to Dirac \cite{dirac_note_1930} and Frenkel \cite{frenkel_wave_1934}.
The main point, in the modern tensor networks formulation, is that instead of solving the Schrödinger equation and then truncating the MPS representation of the quantum state, one can solve the equations of motion projected into a space of restricted bond dimension \cite{haegeman_time-dependent_2011, haegeman_unifying_2016}.

The general formulation of the Dirac-Frenkel Variational Principle~\cite{raab_diracfrenkelmclachlan_2000} is that one looks for a solution $|\varphi\rangle \in \mathcal{M}$ of the Schrödinger equation where $\mathcal{M} \subset \mathcal{H}$ is a manifold of the total Hilbert space $\mathcal{H}$ in which we think that the relevant physical states `live'.

We define $T_{|\varphi\rangle}\mathcal{M}$ the tangent space of $\mathcal{M}$ around the state $|\varphi\rangle$.
The criterion to find $|\varphi\rangle$ is that for every state $|\chi\rangle \in T_{|\varphi\rangle}\mathcal{M}$

```math
    \langle\chi|\left(\frac{\mathrm{d}}{\mathrm{d}t} - \frac{1}{\mathrm{i}\hbar}\hat{H}\right)|\varphi\rangle =0\ ,
```

which can be interpreted as saying that the time evolution procedure should keep $|\varphi\rangle$ inside of the manifold $\mathcal{M}$.

The term *variational* in the name of the method comes from the fact that in practice one aims at minimising the right-hand side of Eq.~(\ref{eq:DiracFrenkel1}) to find $|\varphi\rangle$.

Introducing $\hat{P}_{T_{|\varphi\rangle}\mathcal{M}}$ the projector onto the tangent space $T_{|\varphi\rangle}\mathcal{M}$, we can write the state $|\chi\rangle = \hat{P}_{T_{|\varphi\rangle}\mathcal{M}}|\phi\rangle$ with $|\phi\rangle$ a state in $\mathcal{H}$.
Leading to

```math
    \forall |\phi\rangle \in \mathcal{H}, \ \langle\phi|\hat{P}_{T_{|\varphi\rangle}\mathcal{M}}\left(\frac{\mathrm{d}}{\mathrm{d}t} - \frac{1}{\mathrm{i}\hbar}\hat{H}\right)|\varphi\rangle =0\ .
```
Because the time derivation and the projector commute, we have

```math
    \forall |\phi\rangle \in \mathcal{H}, \ \langle\phi|\left(\frac{\mathrm{d}}{\mathrm{d}t} - \frac{1}{\mathrm{i}\hbar}\hat{P}_{T_{|\varphi\rangle}\mathcal{M}}\hat{H}\right)|\varphi\rangle =0\ .
```
This equation must be true for any $|\phi\rangle \in \mathcal{H}$, Eq.~(\ref{eq:DiracFrenkel1}) can thus be written

```math
    \left(\frac{\mathrm{d}}{\mathrm{d}t} - \frac{1}{\mathrm{i}\hbar}\hat{P}_{T_{|\varphi\rangle}\mathcal{M}}\hat{H}\right)|\varphi\rangle =0\ .
```

In the context of MPS, the manifold $\mathcal{M}$ will correspond to the space of full-ranked MPS of a given bond dimension $D$, and the tangent space will be the space spanned by variations of single MPS tensors.

The major advantage of this method is that it naturally preserves the unitarity of the time evolution and conserves the energy.

## Bibliography
[^chin_exact_2010]:
    > Chin, A. W.; Rivas, Á.; Huelga, S. F.; Plenio, M. B. Exact Mapping between System-Reservoir Quantum Models and Semi-Infinite Discrete Chains Using Orthogonal Polynomials. Journal of Mathematical Physics 2010, 51 (9), 092109. https://doi.org/10.1063/1.3490188.

[^chin_role_2013]:
    > Chin, A. W.; Prior, J.; Rosenbach, R.; Caycedo-Soler, F.; Huelga, S. F.; Plenio, M. B. The Role of Non-Equilibrium Vibrational Structures in Electronic Coherence and Recoherence in Pigment–Protein Complexes. Nature Phys 2013, 9 (2), 113–118. https://doi.org/10.1038/nphys2515.

[^alvertis_nonequilibrium_2019]:
    > Alvertis, A. M.; Schröder, F. A. Y. N.; Chin, A. W. Non-Equilibrium Relaxation of Hot States in Organic Semiconductors: Impact of Mode-Selective Excitation on Charge Transfer. J. Chem. Phys. 2019, 151 (8), 084104. https://doi.org/10.1063/1.5115239.

[^dunnett_influence_2021]:
    > Dunnett, A. J.; Gowland, D.; Isborn, C. M.; Chin, A. W.; Zuehlsdorff, T. J. Influence of Non-Adiabatic Effects on Linear Absorption Spectra in the Condensed Phase: Methylene Blue. J. Chem. Phys. 2021, 155 (14), 144112. https://doi.org/10.1063/5.0062950.

[^caycedosoler_exact_2022]:
    > Caycedo-Soler, F.; Mattioni, A.; Lim, J.; Renger, T.; Huelga, S. F.; Plenio, M. B. Exact Simulation of Pigment-Protein Complexes Unveils Vibronic Renormalization of Electronic Parameters in Ultrafast Spectroscopy. Nat Commun 2022, 13 (1), 2912. https://doi.org/10.1038/s41467-022-30565-4.

[^orus_practical_2014]:
    > Orus, R. A Practical Introduction to Tensor Networks: Matrix Product States and Projected Entangled Pair States. Annals of Physics 2014, 349, 117–158. https://doi.org/10.1016/j.aop.2014.06.013.

[^paeckel_timeevolution_2019]:
    > Paeckel, S.; Köhler, T.; Swoboda, A.; Manmana, S. R.; Schollwöck, U.; Hubig, C. Time-Evolution Methods for Matrix-Product States. Annals of Physics 2019, 411, 167998. https://doi.org/10.1016/j.aop.2019.167998.

[^prior_efficient_2010]:
    > Prior, J.; Chin, A. W.; Huelga, S. F.; Plenio, M. B. Efficient Simulation of Strong System-Environment Interactions. Phys. Rev. Lett. 2010, 105 (5), 050404. https://doi.org/10.1103/PhysRevLett.105.050404.

[^woods_simulating_2015]:
    > Woods, M. P.; Cramer, M.; Plenio, M. B. Simulating Bosonic Baths with Error Bars. Phys. Rev. Lett. 2015, 115 (13), 130401. https://doi.org/10.1103/PhysRevLett.115.130401.

[^tamascelli_efficient_2019]:
    > Tamascelli, D.; Smirne, A.; Lim, J.; Huelga, S. F.; Plenio, M. B. Efficient Simulation of Finite-Temperature Open Quantum Systems. Phys. Rev. Lett. 2019, 123 (9), 090402. https://doi.org/10.1103/PhysRevLett.123.090402.

[^dunnett_simulating_2021]:
    > Dunnett, A. J.; Chin, A. W. Simulating Quantum Vibronic Dynamics at Finite Temperatures With Many Body Wave Functions at 0 K. Front. Chem. 2021, 8. https://doi.org/10.3389/fchem.2020.600731.

[^lacroix_unveiling_2021]:
    > Lacroix, T.; Dunnett, A.; Gribben, D.; Lovett, B. W.; Chin, A. Unveiling Non-Markovian Spacetime Signaling in Open Quantum Systems with Long-Range Tensor Network Dynamics. Phys. Rev. A 2021, 104 (5), 052204. https://doi.org/10.1103/PhysRevA.104.052204.

[^appel_mathematics_2007]:
    > Appel, W. Mathematics for Physics and Physicists; Princeton University Press, 2007.

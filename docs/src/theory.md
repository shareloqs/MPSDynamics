# Theoretical Background

## Chain-Mapping of bosonic environments

We consider, in the Schrödinger picture, a general Hamiltonian where a non-specified system interacts linearly with a bosonic environments

```math
\begin{align}
    \hat{H} =& \hat{H}_S + \int_0^{+\infty} \hbar\omega\hat{a}^\dagger_\omega\hat{a}_\omega\ \mathrm{d}\omega + \hat{A}_S\int_0^{+\infty}\sqrt{J(\omega)}\left(\hat{a}_\omega + \hat{a}^\dagger_\omega\right)\mathrm{d}\omega
\end{align}
```

where ``\hat{a}_\omega`` (``\hat{a}^\dagger_\omega``) is a bosonic annihilation (creation) operator for a normal mode of the environment of energy ``\hbar\omega``, ``\hat{A}_S`` is a system operator, and ``J(\omega) = \sum_k |g_k|^2\delta(\omega - \omega_k)`` is the bath spectral density (SD), defined with the microscopic system-environment coupling strength ``g_k``.
The SD quantifies the coupling strengths of the different normal modes of the environment with the system.
Any SD that is not flat corresponds to a non-Markovian environment.

### Zero Temperature

Let us consider the Hamiltonian presented in Eq.(1).
We can introduce a unitary transformation of the continuous normal modes ``\hat{a}_\omega`` to an infinite discrete set of interacting modes ``\hat{b}_n``[^chin_exact_2010].

```math
\begin{align}
    \hat{a}_\omega &= \sum_{n=0}^{+\infty} U_n(\omega)\hat{b}_n = \sum_{n=0}^{+\infty} \sqrt{J(\omega)}P_n(\omega)\hat{b}_n\ ,
\end{align}
```

where ``P_n(\omega)`` are orthonormal polynomials such that

```math
\begin{align}
    \int_{0}^{+\infty}P_n(\omega)P_m(\omega)J(\omega)\mathrm{d}\omega = \delta_{n,m}\ ;
\end{align}
```

and the inverse transformation is

```math
\begin{align}
    \hat{b}_n &= \int_0^{+\infty} U_n(\omega)\hat{a}_\omega\mathrm{d}\omega\ .
\end{align}
```

Note that the orthonormality of the polynomials ensures the unitarity of the transformation defined in Eq.(2).
The mapping from a continuous set of modes to a (still infinite) discrete set might seem counter-intuitive, however it is a direct consequence of the separability of the underlying Hilbert space.

Under this transformation, the Hamiltonian in Eq.(1) becomes

```math
\begin{align}
    \hat{H}= \hat{H}_S &+ \sum_{n=0}^{+\infty}\varepsilon_n\hat{b}_n^\dagger\hat{b}_n + t_n(\hat{b}_{n+1}^\dagger\hat{b}_n + \mathrm{h.c.}) + \kappa\hat{A}_S(\hat{b}_0 + \hat{b}_0^\dagger)\ .
\end{align}
```

Hence, this mapping transforms the normal bath Hamiltonian into a tight-binding Hamiltonian with on-site energies ``\varepsilon_n`` and hopping energies ``t_n``.
Another important consequence of this mapping is that now the system only interacts with the first mode ``n = 0`` of the chain-mapped environment.
The chain coefficients ``\varepsilon_n``, ``t_n``, and the coupling ``\kappa`` depend solely on the SD.

This makes chain mapping a tool of choice for describing systems coupled to environment with highly structured SD (e.g. experimentally measured or calculated *ab initio*)[^chin_role_2013][^alvertis_nonequilibrium_2019][^dunnett_influence_2021][^caycedosoler_exact_2022]
In this new representation, the Hamiltonian in Eq.(5) has naturally a 1D chain topology.
This makes its representation as a Matrix Product Operator (MPO) and the representation of the joint \{System + Environment\} wave-function as a Matrix Product State (MPS) suited [^orus_practical_2014][^paeckel_timeevolution_2019].

The orthogonal polynomial-based chain mapping and the subsequent representation of the joint wave-function as a MPS (and the operators as MPO) are the building blocks of the Time-dependent Density operator with Orthonormal Polynomials Algorithm (TEDOPA) one of the state-of-the-art numerically exact method to simulate the dynamics of open quantum systems especially in the non-Markovian, non-perturbative regimes both at zero and finite temperatures [^prior_efficient_2010][^woods_simulating_2015][^tamascelli_efficient_2019][^dunnett_simulating_2021][^lacroix_unveiling_2021].

### Finite Temperature

Explain that by extending the bath to negative frequencies and having temperature-dependent system environment couplings, it is possible to describe the finite temperature case as an effective zero temperature one.
Hence, we can keep the pure state description and avoid moving to density matrices at the cost of doubling the size of the environment.

## Computation of the chain coefficient

Explain the three-term recurence relation and the relation of the coefficients to the chain on-site and hopping energies.

Explain that for some weight function/SD they are known analytically and that for others we can use the build-in routines inspired by Gautschi or the PolyChaos.jl package.

## Time-Dependent Variational Principal

Explain the concept behind TDVP and give references

## Tensor Networks

Explain quickly why and give references

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

# Convergence checks

In this section, we analyze the numerical accuracy of our simulations by outlining the necessary approximations and providing practical guidelines to ensure convergence of results.

## Approximations

### Chain length $N$
In numerical simulations, a truncation on the number of chain modes (and therefore the chain length) is introduced to handle a finite chain instead of a semi-infinite one. This truncation, denoted as NN, corresponds to a sampling of the modes in the original environment.

The chain length $N$ is directly connected to the simulation time as follows: excitations injected into the chain propagate as a wavefront traveling along the chain. Perturbations to the initial state outside this wavefront are exponentially suppressed. Therefore, truncating the chain basis beyond the expanding wavefront ensures that the resulting sampling error is also exponentially small  [^woods_simulating_2015][^DeVega_howto_2015].

To optimize the chain length for a given simulation time and set of chain coefficients, the following procedure can be used:

1. Set the total simulation time and compute the chain coefficients for an intentionally oversized chain (longer than needed):
```
T = 100 
N_huge = 1000
cpars = chaincoeffs_flat(N_huge, αchain; ωc = 1.0)
```
2. Use the built-in function `findchainlength` to determine the optimal chain length $N_opt$ based on the propagation speed of the wavefront (given by the hopping coefficients $t_n$):
```
N_opt = findchainlength(T, cpars; eps=10^-4, verbose=false)
```

### Local dimension $d$
Another critical truncation is imposed on the local dimension $d$ of each tensor in the MPS that undergoes dynamic evolution. The local dimension $d$ corresponds to the number of Fock states retained in the Hilbert space of each chain mode. Since harmonic oscillators are, in principle, infinite-dimensional systems, truncating their Hilbert space to a finite $d$ is necessary for numerical computations. The choice of $d$ determines the maximum number of excitations per site in the MPS and must be carefully tuned to capture the relevant physics while ensuring numerical convergence.
For a bath that initially contains only a finite number of particles, i.e. any bath in practice, the error originating from the local Hilbert space truncation vanishes exponentially as $d$ increases[^woods_simulating_2015]. 

A practical guideline for choosing $d$ can be derived by considering the physical states being simulated. For example, for coherent states, the occupation number follows a Poisson distribution. In such cases, the average occupation number is given by the mean $\langle n \rangle$, and the probability of observing a state with occupation number $n$ decreases exponentially for $n \gg \langle n \rangle$. Thus, the truncation $d$ should be chosen such that the cumulative probability of truncation error is negligible.

### Bond dimensions and time evolution
Different time-evolution algorithms are implemented in MPSDynamics, and can be selected as options in the function `runsim`. Detailed descriptions of these algorithms are provided in the *Theoretical Background* section of this documentation. The key parameter for ensuring convergence in these methods is controlled via the `convparams` option, which accepts arrays of multiple parameter values to test until convergence is achieved.

Below are the main time-evolution algorithms and their corresponding convergence parameters:
- One-Site Time-Dependent Variational Principle (1TDVP):
  - Selected using: `method = :1TDVP`.
  - The convergence parameter is the bond dimension of the MPS, which defines the manifold on which the time evolution is constrained. To ensure convergence, multiple bond dimensions can be tested: `convparams = [bond1, bond2, bond3]`.
- Two-Sites Time-Dependent Variational Principle (2TDVP):
  - Selected using: `method = :2TDVP`.
  - The convergence parameter is the SVD singular value truncation threshold, which determines the precision of the truncation during simulations. Multiple truncation thresholds can be tested: `convparams = [trunc1, trunc2, trunc3]`.
  - The bond dimensions at each chain-site and at each time-step can be saved by specifying `savebonddims = true`. The maximum allowed value of bond dimension can be chosen by setting the `Dlim` option.
- Adaptive Time-Dependent Variational Principle (DTDVP):
  - Selected using: `method = :DTDVP`.
  - The convergence parameter is the target precision for the adaptive algorithm. Multiple precision values can be provided to test convergence: `convparams = [prec1, prec2, prec3]`.
  - The bond dimensions at each chain-site and at each time-step can be saved by specifying `savebonddims = true`. The maximum allowed value of bond dimension can be chosen with the `Dlim` option.
  - In some cases, it might be better to embed the initial time MPS in a manifold of higher bond-dimension, which can be done with the function `mpsembed!`, to avoid for the dynamics to get stuck.
 
### Time-step 
The time step used in simulations must be small enough to accurately capture the dynamics of the system. To ensure convergence, gradually reduce the time step until key observables (e.g., energy, population dynamics, or correlation functions) stabilize. It is often practical to test a range of time steps and assess their impact on results to determine an optimal balance between accuracy and computational cost.

## Common pitfalls

### Hard cut-off in the spectral density function
The cutoff frequency $\omega_c$ in the spectral density function has to be selected in order to incorporate all of the relevant frequencies. It is to be noted that, for spectral densities belonging to the Szego class[^chin_exact_2010], the cutoff frenquency $\omega_c$ also determines the asymptotic values to which the chain coefficients converge [^woods_mappings_2014]


### MPS gauge choice 

# References
[^woods_simulating_2015]: Woods, M. P.; Cramer, M.; Plenio, M. B. Simulating Bosonic Baths with Error Bars. Phys. Rev. Lett. 2015, 115 (13), 130401. https://doi.org/10.1103/PhysRevLett.115.130401.
[^DeVega_howto_2015]: De Vega, I.; Schollwöck, U.; Wolf, F. A. How to Discretize a Quantum Bath for Real-Time Evolution. Phys. Rev. B 2015, 92 (15), 155126. https://doi.org/10.1103/PhysRevB.92.155126.
[^woods_mappings_2014]: Woods,  M. P. et al. “Mappings of open quantum systems onto chain representations and Markovian embeddings”. In: Journal of Mathematical Physics 55.3 (Mar. 2014), p. 032101. issn: 0022-2488, 1089-7658. doi: 10.1063/1.4866769
[^chin_exact_2010]: Chin, A. W.  et al. “Exact mapping between system-reservoir quantum models and semi-infinite discrete chains using orthogonal polynomials”. In: Journal of Mathematical Physics 51.9 (Sept. 2010), p. 092109. issn: 0022-2488, 1089-7658. doi: 10.1063/1.3490188. arXiv: 1006.4507

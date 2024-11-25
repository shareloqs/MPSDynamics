# Convergence checks

In this section, we analyze the numerical accuracy of our simulations by outlining the necessary approximations and providing practical guidelines to ensure convergence of results.

## Approximations

### Chain length $N$
In numerical simulations, a truncation on the number of chain modes (and therefore the chain length) is introduced to handle a finite chain instead of a semi-infinite one. This truncation, denoted as NN, corresponds to a sampling of the modes in the original environment.

The chain length $N$ is directly connected to the simulation time as follows: excitations injected into the chain propagate as a wavefront traveling along the chain. Perturbations to the initial state outside this wavefront are exponentially suppressed. Therefore, truncating the chain basis beyond the expanding wavefront ensures that the resulting sampling error is also exponentially small.

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
Another critical truncation is imposed on the local dimension $d$ of each tensor in the MPS that undergoes dynamic evolution. The local dimension dd corresponds to the number of Fock states retained in the Hilbert space of each chain mode. Since harmonic oscillators are, in principle, infinite-dimensional systems, truncating their Hilbert space to a finite $d$ is necessary for numerical computations. The choice of $d$ determines the maximum number of excitations per site in the MPS and must be carefully tuned to capture the relevant physics while ensuring numerical convergence.

A practical guideline for choosing $d$ can be derived by considering the physical states being simulated. For example, for coherent states, the occupation number follows a Poisson distribution. In such cases, the average occupation number is given by the mean $\langle n \rangle$, and the probability of observing a state with occupation number $n$ decreases exponentially for $n \gg \langle n \rangle$. Thus, the truncation $d$ should be chosen such that the cumulative probability of truncation error is negligible.

### Time-step error
Different time-evolution algorithms are implemented in MPSDynamics, and can be selected as options in the function `runsim`. Detailed descriptions of these algorithms are provided in the *Theoretical Background* section of the documentation. The key parameter for ensuring convergence in these methods is controlled via the `convparams` option, which accepts arrays of multiple parameter values to test until convergence is achieved.

Below are the main time-evolution algorithms and their corresponding convergence parameters:
- One-Site Time-Dependent Variational Principle (1TDVP):
  - Selected using: `method = :1TDVP`.
  - The convergence parameter is the bond dimension of the MPS, which defines the manifold on which the time evolution is constrained. To ensure convergence, multiple bond dimensions can be tested: `convparams = [bond1, bond2, bond3]`.
- Two-Sites Time-Dependent Variational Principle (2TDVP):
  - Selected using: `method = :2TDVP`.
  - The convergence parameter is the SVD singular value truncation threshold, which determines the precision of the truncation during simulations. Multiple truncation thresholds can be tested: `convparams = [trunc1, trunc2, trunc3]`
- Adaptive Time-Dependent Variational Principle (DTDVP):
  - Selected using: `method = :DTDVP`.
  - The convergence parameter is the target precision for the adaptive algorithm. Multiple precision values can be provided to test convergence: `convparams = [prec1, prec2, prec3]`

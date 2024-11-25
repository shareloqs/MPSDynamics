# Convergence checks

In this section, we analyze the numerical accuracy of our simulations by outlining the necessary approximations and providing practical guidelines to ensure convergence of results.

### Chain length $N$
During a numerical simulation, a truncation on the number of chain modes (and therefore chain length) will be introduced, in order to work with a chain of finite length instead of a semi-infinite one. This truncation on chain modes, let us say $N$, corresponds to a sampling on the modes in the original environment. 

The length chain is connected to the simulation time in the following way: excitations injected from the system form a wave-front that travels along the chain. Any perturbation to the initial state outside of this wave-front are exponentially suppressed. It is thus natural to truncate in the chain basis, making sure that this truncation happens **beyond** the expanding wave-front, as in this way the sampling error is also exponentially small [^woods_simulating_2015][^DeVega_howto_2015]. 

To make sure that, for a given simulation time, and for a given set of chain coefficients, the chain-length is the optimal one, we can do as follows. We set the total simulation time, and we compute the chain coefficients for a non-optimal chain, much longer than needed.
```
T = 100 
N_huge = 1000
cpars = chaincoeffs_flat(N_huge, αchain; ωc = 1.0)
```
We can then use the built-in function `findchainlength`, that computes the required chain length `N_opt` from the propagation speed on the chain sites (given by the hopping coefficients $t_n$):
```
N_opt = findchainlength(T, cpars; eps=10^-4, verbose=false)
```

### Local dimension $d$
Another critical truncation is imposed on the local dimension $d$ of each tensor in the MPS that undergoes dynamic evolution. The local dimension $d$ corresponds to the number of Fock states retained in the Hilbert space of each chain mode. Since harmonic oscillators are, in principle, infinite-dimensional systems, truncating their Hilbert space to a finite $d$ is necessary for numerical computations. The choice of $d$ determines the maximum number of excitations per site in the MPS and must be carefully tuned to capture the relevant physics while ensuring numerical convergence.
For a bath that initially contains only a finite number of particles, i.e. any bath in practice, the error originating from the local Hilbert space truncation vanishes exponentially as $d$ increases[^woods_simulating_2015]. 

### 


# References
[^woods_simulating_2015]: Woods, M. P.; Cramer, M.; Plenio, M. B. Simulating Bosonic Baths with Error Bars. Phys. Rev. Lett. 2015, 115 (13), 130401. https://doi.org/10.1103/PhysRevLett.115.130401.
[^DeVega_howto_2015]: De Vega, I.; Schollwöck, U.; Wolf, F. A. How to Discretize a Quantum Bath for Real-Time Evolution. Phys. Rev. B 2015, 92 (15), 155126. https://doi.org/10.1103/PhysRevB.92.155126.

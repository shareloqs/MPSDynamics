# Convergence checks

In this section, we analyze the numerical accuracy of our simulations by outlining the necessary approximations and providing practical guidelines to ensure convergence of results.

## Chain length $N$
During a numerical simulation, a truncation on the number of chain modes (and therefore chain length) will be introduced, in order to work with a chain of finite length instead of a semi-infinite one. This truncation on chain modes, let us say $N$, introduces a sampling on the modes in the original environment. 

## Local dimension $d$
Another critical truncation is imposed on the local dimension $d$ of each tensor in the MPS that undergoes dynamic evolution. The local dimension dd corresponds to the number of Fock states retained in the Hilbert space of each chain mode. Since harmonic oscillators are, in principle, infinite-dimensional systems, truncating their Hilbert space to a finite $d$ is necessary for numerical computations. The choice of dd determines the maximum number of excitations per site in the MPS and must be carefully tuned to capture the relevant physics while ensuring numerical convergence.

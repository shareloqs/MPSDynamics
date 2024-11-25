# Convergence checks

In this page we give some tips on how to obtain converged simulations. 


## Chain length $N$
During a numerical simulation, a truncation on the number of chain modes (and therefore chain length) will be introduced, in order to work with a chain of finite length instead of a semi-infinite one. This truncation on chain modes, let us say $N$, introduces a sampling on the modes in the original environment.

## Local dimension $d$
A second truncation has to be introduced on the local dimension $d$ of each tensor comprising the MPS that will be dynamically evolved. The local dimension coincides with the truncated number of Fock's states of the chain mode's Hilbert space: each site of the MPS, can be filled up by at most $d$ excitations. Since in this case the excitations are bosonic, the Pauli exclusion principle does not hold, and this number should in principle be infinite.

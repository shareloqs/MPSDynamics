# MPSDynamics.jl

Tensor network simulations for finite temperature, open quantum system dynamics.

This package is intended to provide an easy to use interface for performing tensor network simulations on Matrix Product
States (MPS). MPSDynamics.jl is a versatile package which supports both chain and (loop-free) tree MPS, as well as
providing a choice of several time evolution algorithms. The package also provides strong support for the measurement
of observables, as well as the storing and logging of data, which makes it a useful tool for the study of many-body
physics. The package was originally developped with the aim of studying open system dynamics at finite temperature using
the T-TEDOPA mapping [1], however the methods implemented can equally be applied to other areas of physics.

The methods currently implemented are

* 1-site TDVP on tree and chain MPS [2]
* 2-site TDVP on chain MPS [2]
* a variant of 1-site TDVP with dynamic bond-dimensions on chain MPS [3]

The elementary tensor operations are implemented in all cases using the [TensorOperations.jl](https://github.com/Jutho/TensorOperations.jl) package.

# Installation

The package may be installed by typing the following into a Julia REPL

```julia
] add https://github.com/angusdunnett/MPSDynamics.git
```

# Usage

The basic usage is as follows. First, include the package.

```julia
using MPSDynamics
```

To set up a simulation we require an MPS representing our initial wave-function and a Matrix Product Operator (MPO) representing our Hamiltonian.

MPSDynamics.jl contains various functions for generating MPSs and MPOs used for simulating certain models, but no attempt is made to be comprehensive. For generic MPO constuction, one can use the [ITensors.jl](https://github.com/ITensor/ITensors.jl) package and convert the resulting object into a form compatible with MPSDynamics.jl using the function MPOtoVector.

In this example 



# Publications

* Simulating quantum vibronic dynamics at finite temperatures with many body wave functions at 0K
     * [https://doi.org/10.3389/fchem.2020.600731](https://doi.org/10.3389/fchem.2020.600731)

* Matrix Product State Simulations of Non-Equilibrium Steady States and Transient Heat Flows in the Two-Bath Spin-Boson Model at Finite Temperatures
     * [https://doi.org/10.3390/e23010077](https://doi.org/10.3390/e23010077)

# Data Repositories

* Exact Spin-Boson-Model Tunnelling Dynamics with Time Dependent Variation Matrix Product States (TDVMPS). Barrier height and temperature parameter space.
     * [10.5281/zenodo.4432014](https://doi.org/10.5281/zenodo.4432014)

* Real-time benchmark dynamics of the Ohmic Spin-Boson Model computed with Time-Dependent Variational Matrix Product States. (TDVMPS) coupling strength and temperature parameter space.
     * [10.5281/zenodo.4352728](https://doi.org/10.5281/zenodo.4352728)

# References

* [1] [D. Tamascelli, A. Smirne, J. Lim, S. F. Huegla, and M. B. Plenio, Physical Review Letters 123, 090402 (2019) arXiv: 1811.12418](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.123.090402)

* [2] [J. Haegeman, C. Lubich, I. Oseledets, B. Vandereycken, and F. Verstraete, Physical Review B 94, 165116 (2016), arXiv: 1408.5056](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.94.165116)

* [3] [A. J. Dunnett, and A. W. Chin, arXiv : 2007.13528](https://arxiv.org/abs/2007.13528)
# MPSDynamics.jl

Tensor network simulations for finite temperature, open quantum system dynamics.

This package is intended to provide an easy to use interface for performing tensor network simulations on Matrix Product
States (MPS). MPSDynamics.jl is a versatile package which supports both chain and (loop-free) tree MPS, as well as
providing a choice of several time evolution algorithms. The package also provides strong support for the measurement
of observables, as well as the storing and logging of data, which makes it a useful tool for the study of many-body
physics. The package was originally developed with the aim of studying open system dynamics at finite temperature using
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

MPSDynamics.jl contains various functions for generating MPSs and MPOs used for simulating certain models, but no attempt is made to be comprehensive. For generic MPO construction, one can use the [ITensors.jl](https://github.com/ITensor/ITensors.jl) package and convert the resulting object into a form compatible with MPSDynamics.jl using the function `MPOtoVector`.

In this example we will consider the spin-boson model. First we define parameters and generate the MPO.

```julia
d=6
N=30

α = 0.5
Δ = 0.0
ω0 = 0.2
s = 1
cpars = chaincoeffs_ohmic(N, α, s)

H = spinbosonmpo(ω0, Δ, d, N, cpars)

```

Then we create the MPS.

```julia
A = productstatemps(physdims(H));
```

This will generate a product state MPS with local Hilbert space dimensions corresponding to the MPO `H`, representing the
spin in the up state and the bath in the vacuum state.

We may then wish to construct some observables to measure along the trajectory. For example

```julia
ob1 = OneSiteObservable("sz", sz, 1)
```

creates an object which represents the measurement of the expectation <img
src="https://render.githubusercontent.com/render/math?math=\sigma_x"> on the first site of the chain, i.e. on the spin. The
string passed to the first argument is just a label that will be used to retrieve the measurement data after the run.

We may also wish to measure the bath observables.

```julia
ob2 = OneSiteObservable("chain mode occupation", numb(d), (2,N+1))
```

This will measure the number operator (truncated to d Fock states) on all chain modes, i.e. on sites 2 to
N+1 inclusive.

Finally we launch the simulation with the function `runsim`.

```julia
dt = 0.2
T = 60.0

A, dat = runsim(dt, T, A, H;
                name = "ohmic spin boson model",
                method = :TDVP1,
                obs = [ob2],
                convobs = [ob1],
                params = @LogParams(N, d, α, Δ, ω0, s),
                convparams = [2,4,8],
                verbose = false,
                save = true,
                plot = true,
                );
```

This will propagate the MPS up to time `T` in time steps of `dt`. The simulation will be performed using 1-site TDVP with
bond-dimensions of 2, 4 and 6 in order to check for convergence. The observables supplied to `convobs` will be measured
at every time step for every bond-dimension, while the observables supplied to `obs` will only be measured for the last
(most accurate) convergence parameter supplied to `convparams`.

The final MPS is returned to `A` and the measurement data is returned to `dat`. If the option `save=true` is used the
data will also be saved to a file. The save directory may be specified using the option `savedir`, by default the save
directory is ~/MPSDynamics, which will be created if it doesn't exist (if using Windows the slashes will need to be reversed).

The data is stored in the JLD format which is based on HDF5. Loading the data in julia using the
[JLD](https://github.com/JuliaIO/JLD.jl) package will recover the full type information of the Julia variables that were
stored. At the same time the HDF5 format is natively supported across many platforms and languages. For example, to load
the spin-z data in Mathematica one could do

```mathematica
Import["~/MPSDynamics/XXXXX/dat_XXXXX.jld",{"HDF5","Datasets","/data/sz"}]
```



# Publications

* Simulating quantum vibronic dynamics at finite temperatures with many body wave functions at 0K
     * [https://doi.org/10.3389/fchem.2020.600731](https://doi.org/10.3389/fchem.2020.600731)

* Matrix Product State Simulations of Non-Equilibrium Steady States and Transient Heat Flows in the Two-Bath Spin-Boson Model at Finite Temperatures
     * [https://doi.org/10.3390/e23010077](https://doi.org/10.3390/e23010077)

# Data Repositories

* Exact Spin-Boson-Model Tunneling Dynamics with Time Dependent Variation Matrix Product States (TDVMPS). Barrier height and temperature parameter space.
     * [10.5281/zenodo.4432014](https://doi.org/10.5281/zenodo.4432014)

* Real-time benchmark dynamics of the Ohmic Spin-Boson Model computed with Time-Dependent Variational Matrix Product States. (TDVMPS) coupling strength and temperature parameter space.
     * [10.5281/zenodo.4352728](https://doi.org/10.5281/zenodo.4352728)

# References

* [[1]](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.123.090402) D. Tamascelli, A. Smirne, J. Lim, S. F. Huegla, and M. B. Plenio, Physical Review Letters 123, 090402 (2019) arXiv: 1811.12418

* [[2]](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.94.165116) J. Haegeman, C. Lubich, I. Oseledets, B. Vandereycken, and F. Verstraete, Physical Review B 94, 165116 (2016), arXiv: 1408.5056

* [[3]](https://arxiv.org/abs/2007.13528) A. J. Dunnett, and A. W. Chin, arXiv : 2007.13528
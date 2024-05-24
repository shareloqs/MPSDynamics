# Introduction

The `MPSDynamics.jl` package provides an easy to use interface for performing tensor network simulations on matrix product states (MPS) and tree tensor network (TTN) states.
Written in the Julia programming language, `MPSDynamics.jl` is a versatile open-source package providing a choice of several variants of the Time-Dependent Variational Principle (TDVP) method for time evolution. 
The package also provides strong support for the measurement of observables, as well as the storing and logging of data, which makes it a useful tool for the study of many-body physics. 
The package has been developed with the aim of studying non-Markovian open system dynamics at finite temperature using the state-of-the-art numerically exact Thermalized-Time Evolving Density operator with Orthonormal Polynomials Algorithm (T-TEDOPA) based on environment chain mapping.
However the methods implemented can equally be applied to other areas of physics.

!!! warning
	The documentation is currently undergoing massive restructurations/improvement. It's a work in progress until the next release scheduled for May, 2024.

## Installation

The package may be installed by typing the following into a Julia REPL

```julia
    ] add https://github.com/shareloqs/MPSDynamics.git
```

## Table of Contents

```@contents
Pages = ["index.md", "user-guide.md", "examples/sbm.md", "examples/puredephasing.md", "examples/timedep.md", "examples/anderson-model.md", "examples/bath-observables.md", "examples/protontransfer.md", "theory.md", "methods.md", "dev.md"]
Depth = 3
```

## Citation
If you use the package in your research, please consider citing it.
You can add the Zenodo record to your BibTex file:

```tex
@misc{mpsdynamics_zenodo2021,
	title = {shareloqs/{MPSDynamics}},
	shorttitle = {{MPSDynamics.jl}},
	url = {https://zenodo.org/record/5106435},
	abstract = {Tensor network simulations for finite temperature, open quantum system dynamics},
	publisher = {Zenodo},
	author = {Dunnett, Angus and Lacroix, Thibaut and Le Dé, Brieuc and Riva, Angela},
	year = {2021},
	doi = {10.5281/zenodo.5106435},
}
```

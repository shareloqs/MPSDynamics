# Introduction

The `MPSDynamics.jl` package provides an easy to use interface for performing tensor network simulations on matrix product states (MPS) and tree tensor network (TTN) states.
Written in the Julia programming language, `MPSDynamics.jl` is a versatile open-source package providing a choice of several variants of the Time-Dependent Variational Principle (TDVP) method for time evolution. 
The package also provides strong support for the measurement of observables, as well as the storing and logging of data, which makes it a useful tool for the study of many-body physics. 
The package has been developed with the aim of studying non-Markovian open system dynamics at finite temperature using the state-of-the-art numerically exact Thermalized-Time Evolving Density operator with Orthonormal Polynomials Algorithm (T-TEDOPA) based on environment chain mapping.
However the methods implemented can equally be applied to other areas of physics.

## Table of Contents

```@contents
```

## Installation

The package may be installed by typing the following into a Julia REPL

```julia
    ] add https://github.com/shareloqs/MPSDynamics.git
```

## Functions

```@autodocs
Modules = [MPSDynamics]
```

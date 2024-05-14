# Inspecting the bath by undoing the chain mapping

Here we give some context on the example script provided in `MPSDynamics/example/bath-observables.jl`. This example demonstrates the setup and execution of a simulation for a two-level system coupled to an Ohmic bath at finite temperature, where we exploit the access to the chain observables to:
- undo the chain mapping, thus obtaining their representation in the extended bath of T-TEDOPA, characterized by $J(\omega, \beta)$
- inverting the thermofield transformation, thus obtaining the representation of the physical frequencies in the original environment, characterized by $J(\omega)$

T-TEDOPA allows to substitute a thermally occupied bath by an extended one, in the pure state of the vacuum, extending the bath of frequencies to negative values: the _creation_ from the system of a mode of _negative frequency_ in the extended bath of frequencies corresponds to the _absorption_ of energy for the system from the thermally occupied modes of the environment. 
Recovering the picture of the occupations of the steady states in the physical bath of frequencies is therefore an important step to complete the description.

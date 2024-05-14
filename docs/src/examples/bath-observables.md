# Inspecting the bath by undoing the chain mapping

Here we give some context on the example script provided in `MPSDynamics/example/bath-observables.jl`. This example demonstrates the setup and execution of a simulation for a two-level system coupled to an Ohmic bath at zero temperature, where we exploit the access to the chain observables to:
- undo the chain mapping, thus obtaining their representation in the extended bath of T-TEDOPA, characterized by $J(\omega, \beta)$
- inverting the thermofield transformation, thus obtaining the representation of the physical frequencies in the original environment, characterized by $J(\omega)$

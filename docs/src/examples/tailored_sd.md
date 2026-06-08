# Tailoring the Spectral Density

## Context

Here we detail an example to tailor the Spectral Density (SD) used to represent the bath. While the initial chain mapping procedure has been automatized for thermalized Ohmic type of SD (thanks to the methods [`MPSDynamics.chaincoeffs_ohmic`](@ref) and [`MPSDynamics.chaincoeffs_finiteT`](@ref)), it is possible to describe other types of SD. In this example, the SD is treated either with a continuous definition for several frequency domains or with a discrete description with frequency values and corresponding spectral density values as inputs. 

The Spin-Boson Model is chosen to illustrate this example. For details about the model, we refer the interested reader to the example [The Spin-Boson Model](@ref).

## The code to define a specific spectral density

We describe here only differences with the Spin-Boson Model in order to parametrize the Hamiltonian with a specific spectral density.

First, we load the `MPSdynamics.jl` package to be able to perform the simulation, the `Plots.jl` one to plot the results, the `LaTeXStrings.jl` one to be able to use ``\LaTeX`` in the plots. A notation is introduced for `+∞`.

```julia
using MPSDynamics, Plots, LaTeXStrings

const ∞  = Inf
```
We then define variables for the physical parameters of the simulation.

```julia
#----------------------------
# Physical parameters
#----------------------------

ω0 = 0.2 # TLS gap

Δ = 0.0 # tunneling 
``` 

Now comes the definition of the spectral density. The type of SD has to be chosen between `:Discrete` and `:Continuous` (one choice has to be commented instead of the other). 

For the `:Discrete` definition, the bath inputs become the discrete frequencies (`freqs`) with the corresponding spectral density values (`Jω`). 

For the `:Continuous` definition, the bath inputs are as usual. Among these, two are convergence parameters:

* `d` is the number of states we retain for the truncated harmonic oscillators representation of environmental modes
* `N` is the number of chain (environmental) modes we keep. This parameters determines the maximum simulation time of the simulation: indeed excitations that arrive at the end of the chain are reflected towards the system and can lead to unphysical results

The value of `β` determines whether the environment is thermalized or not. The example as it is is the zero-temperature case. For the finite-temperature case, 'β = ∞' has be commented instead of the line above that precises a `β` value.

After the parameters, the SD formulation has still to be defined according to the parameters. For that purpose, frequency domains have to be defined with the `AB` list and a function has to describe the spectral density formula for each domain. In order to follow the T-TEDOPA formulation, formula have to be modified as illustrated when the environment is thermalized. The presented example reproduces an Ohmic SD but frequency domains and formula can be changed at will. 

```julia
#-----------------------
# Options and parameters for the Spectral Density
# Definition of the frequency ranges and J(ω) formula for the continuous case
#-----------------------

#Jω_type=:Discrete

Jω_type=:Continuous

if Jω_type==:Discrete
    
    freqs = [0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50] # Frequency value
    
    Jω = [0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.10] # Value of the spectral density at the respective frequency
    
    N = length(freqs) # The number of mode in the chain is set with the number of discrete frequencies
    
    d = 6 # number of Fock states of the chain modes

elseif Jω_type==:Continuous 

    N = 30 # length of the chain

    d = 6 # number of Fock states of the chain modes

    α = 0.1 # coupling strength
    
    s = 1 # ohmicity
    
    ωc = 1.0 # Cut-off of the spectral density J(ω) 
    
    #β = 100 # Thermalized environment
    β = ∞ # Case zero temperature T=0, β → ∞

    # AB segments correspond to the different frequency ranges where the spectral density is defined.
    # It corresponds to i=1 ; i=2 ; i=3 and i=4 of the Jω_fct(ω,i) function.

    # The y formula and the frequency ranges can be changed while keeping the extra terms for thermalized environments.
    # This example reproduces an Ohmic type spectral density. 

    AB = [[-Inf -ωc];[-ωc 0];[0 ωc];[ωc Inf]]

    function Jω_fct(ω,i)
        if i==1 # First segment of AB (here [-Inf -ωc])
            y = 0
        elseif i==2 # Second segment of AB (here [-ωc 0])
            if β == ∞
                y = 0 # zero for negative frequencies when β == Inf   
            else
		y = -2*α*abs.(ω).^s / ωc^(s-1) .* (coth.((β/2).*ω) .+ 1) ./2 # .* (coth.((β/2).*ω) .+ 1) ./2 and the abs.(ω) have to be added when β != Inf (T-TEDOPA formulation) 
            end
        elseif i==3 # Third segment of AB (here [0 ωc])
            if β == ∞
                y = 2*α*ω.^s / ωc^(s-1) 
            else
		y = 2*α*ω.^s / ωc^(s-1) .* (coth.((β/2).*ω) .+ 1) ./2  # .* (coth.((β/2).*ω) .+ 1) ./2 has to be added when β != Inf (T-TEDOPA formulation) 
            end       
        elseif i==4 # Fourth segment of AB (here [ωc Inf])
            y = 0
        end
        return y
    end
    
else

    throw(ErrorException("The spectral density has either to be Continuous or Discrete"))

end
```

Now that the type of SD and the corresponding parameters have been filled in, the chain parameters have to be calculated. 

For both cases, the method [`MPSDynamics.chaincoeffs_finiteT`](@ref) can be called with the respective inputs.

```julia
#-----------------------
# Calculation of chain parameters
#-----------------------

if Jω_type==:Discrete

    cpars = chaincoeffs_finiteT_discrete(β, freqs, Jω; procedure=:Lanczos, Mmax=5000, save=false)  # chain parameters, i.e. on-site energies ϵ_i, hopping energies t_i, and system-chain coupling c_0

#= If cpars is stored in "../ChainOhmT/discreteT" 
    curdir = @__DIR__
    dir_chaincoeff = abspath(joinpath(curdir, "../ChainOhmT/discreteT"))
    cpars  = readchaincoeffs("$dir_chaincoeff/chaincoeffs.h5", N, β) # chain parameters, i.e. on-site energies ϵ_i, hopping energies t_i, and system-chain coupling c_0
=#

elseif Jω_type==:Continuous

    cpars = chaincoeffs_finiteT(N, β, false; α=α, s=s, J=Jω_fct, ωc=ωc, mc=size(AB)[1], mp=0, AB=AB, iq=1, idelta=2, procedure=:Lanczos, Mmax=5000, save=false)  # chain parameters, i.e. on-site energies ϵ_i, hopping energies t_i, and system-chain coupling c_0

end
```

Next, similarly to the spin-boson example, the simulation parameters are set with the initial MPO and MPS in order to carry out the dynamics propagation. Therefore, this modified example shows an easy way to tailor the spectral density for any system and Hamiltonian.

## Dealing with diverging spectral densities

Here, we explain how the `jacobi` parameter of [`MPSDynamics.chaincoeffs_finiteT`](@ref) can help in the case of diverging (but integrable) spectral densities. The most commonly encountered case is for sub-Ohmic (``s<1``) spectral densities at finite temperature, where due to the ``\coth(\beta\omega/2)`` term, the temperature-dependent spectral density scales a ``J_\beta(\omega) \sim |\omega|^{s-1}`` near ``\omega=0``. The `jacobi` parameter allows you to use a specific integrator that performs better for these types of divergences.

The `jacobi` parameter takes in a matrix with the same shape as `AB`. Each row of `jacobi` corresponds to an interval of `AB`, and the first and second column correspond to the scaling exponent of the spectral density near the left and right interval boundary respectively (should be 0 otherwise). This is only available for finite intervals, so it might be necessary to further split the spectral density support in more `AB` intervals.

Taking the same spectral density as the continuous case of [The code to define a specific spectral density](@ref), it scales as ``\sim \omega^{s-1}`` both near the right boundary of the interval `AB[2]=[-ωc 0]` and near the left boundary of interval `AB[3]=[0 ωc]`. Thus the code to generate the chain coefficients is slightly modified as such.

```julia
#-----------------------
# Calculation of chain parameters
#-----------------------

@assert Jω_type==:Continuous "use the continuous spectral density"
@assert !isinf(β) "should be the finite-temperature case (use a finite β)"

jacobi = [[0 0];[0 s-1];[s-1 0];[0 0]]

cpars = chaincoeffs_finiteT(N, β, false; α=α, s=s, J=Jω_fct, ωc=ωc, mc=size(AB)[1], mp=0, AB=AB, iq=1, idelta=2, procedure=:Lanczos, Mmax=5000, save=false, jacobi=jacobi)  # chain parameters, i.e. on-site energies ϵ_i, hopping energies t_i, and system-chain coupling c_0
```

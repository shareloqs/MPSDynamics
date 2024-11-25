# List of all methods

## Matrix Product State (MPS)

```@autodocs
Modules = [MPSDynamics]
Pages = ["mpsBasics.jl", "reshape.jl","tensorOps.jl", "fundamentals.jl", "switchmpo"]
```
## Tree Tensor Network (TTN)

```@autodocs
Modules = [MPSDynamics]
Pages = ["treeBasics.jl", "treeTDVP.jl","treeMeasure.jl", "treeIterators.jl", "treeTDVP.jl"]
```

## Measure and Obervables 

```@autodocs
Modules = [MPSDynamics]
Pages = ["measure.jl", "observables.jl"]
```

## Models and Hamiltonians (MPO) 

```@autodocs
Modules = [MPSDynamics]
Pages = ["models.jl", "electronkmpo.jl", "excitonphononmpo.jl"]
```

## Chain-Mapping

```@autodocs
Modules = [MPSDynamics]
Pages = ["finitetemperature.jl", "chainOhmT.jl", "gauss.jl", "lanczos.jl", "mcdis2.jl" , "ohmicSpec.jl" , "quadohmT.jl", "r_jacobi.jl", "stieltjes.jl", ]
```

## Dynamics propagation function

```@autodocs
Modules = [MPSDynamics]
Pages = ["MPSDynamics.jl", "chain2TDVP.jl", "chainA1TDVP.jl","chainDMRG.jl","chainTDVP.jl", "run_1TDVP.jl", "run_1TDVPLC.jl", "run_2TDVP.jl", "run_all.jl","run_DTDVP.jl"]
```

## Advanced
```@autodocs
Modules = [MPSDynamics]
Pages = ["flattendict.jl", "logging.jl", "logiter.jl","machines.jl"]
```

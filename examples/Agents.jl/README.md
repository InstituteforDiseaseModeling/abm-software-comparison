# Agents.jl
Agents.jl is a pure Julia framework for agent-based modeling (ABM): a computational simulation methodology where autonomous agents react to their environment (including other agents) given a predefined set of rules. 

## Setup
After installing julia, navigate to the `examples/Agents.jl` directory:
```
julia
using Pkg
Pkg.activate()
```
if this is the first time using the package you'll need to install the dependencies:
```
Pkg.instantiate()
```
otherwise, you can skip to importing a model e.g.,:
```
include("ModelA.jl")
import .ModelA
ModelA.run_model()
```

- [git repository](https://github.com/JuliaDynamics/Agents.jl)
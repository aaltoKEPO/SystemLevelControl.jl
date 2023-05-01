```@meta
CurrentModule = SystemLevelControl
```

```@raw html
<p align="center"><img alt="Logo" src="assets/logo_SLS_bg.png"></p>
<h4 align="center">Welcome to the documentation for SystemLevelControl.jl!</h4>
```

---

## Introduction

__SystemLevelControl.jl__ is a Julia toolbox for synthesizing controllers using the System Level Synthesis (SLS) methodology. The package provides a straightforward, `@Distributed`-enabled, interface for optimal and robust control of large-scale cyberphysical systems. 

We are still not on the official registry. To install __SystemLevelControl.jl__, run 
```julia
import Pkg; Pkg.add(url="https://github.com/aaltoKEPO/SystemLevelControl.jl")
```

The documentation is organized in three sections:
- The [Manual](manual/index.md) details the main concepts and methods implemented in the package. 
- The [Examples](examples/index.md) section provides some tutorials on using the package for solving relevant benchmark problems.
- The [Reference](reference/index.md) section contains the documentation of all important types and functions from the library.

## Notes

!!! note 
    This package is currently under active development. A stable version is to be released soon.

Although general-purpose, the SLS methodology was designed for large-scale cyberphysical systems with sparse communication and actuation networks. This is reflected on our design choices for the package: __SystemLevelControl.jl__ focus on linear discrete-time systems represented through `SparseArrays` data-types. The linear fractional transformation (LFT) framework is used as the theoretical backbone of all methods.

We hope to improve __SystemLevelControl.jl__ until it can serve as a general-purpose control framework (while still bringing all the power of SLS). In the meanwhile, we can recommend other excellent Julia packages for doing analysis or solving more general control problems:

- [`ControlSystems.jl`](https://juliacontrol.github.io/ControlSystems.jl/stable/) and its associated Ecosystem provide a wide collection of tools for analysis and design of control systems. It provides a similar interface as that of the popular Control Systems Toolbox in MATLAB®.
- [`JuMP.jl`](https://jump.dev/JuMP.jl/stable/) is perhaps the most popular modelling framework for mathematical optimization in Julia. A wide class of optimal control problems (including LMIs/SDPs) can be solved with this package. **SystemLevelControl.jl** actually use JuMP to solve SLS problems for the general cases.
- [`TrajectoryOptimization.jl`](http://roboticexplorationlab.org/TrajectoryOptimization.jl/stable/) is a popular framework for solving trajectory optimization problems in Julia, specially for applications in robotics. A distinct feature is the possibility to easily model nonlinear control problems.
- [`JuliaSimControl.jl`](https://help.juliahub.com/juliasimcontrol/stable/) is the package within the JuliaSim Ecosystem that allows for the modelling, analysis and deployment of control systems in a centralized package.

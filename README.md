# DftFunctionals

| **Documentation** | **Build status** | **License** |
|:----------------- |:---------------- |:----------- |
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://juliamolsim.github.io/DftFunctionals.jl/dev) [![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://juliamolsim.github.io/DftFunctionals.jl/stable) | [![Build Status](https://github.com/JuliaMolSim/DftFunctionals.jl/actions/workflows/CI.yml/badge.svg?branch=master)](https://github.com/JuliaMolSim/DftFunctionals.jl/actions/workflows/CI.yml?query=branch%3Amaster) [![Coverage](https://codecov.io/gh/JuliaMolSim/DftFunctionals.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/JuliaMolSim/DftFunctionals.jl) | [![][license-img]][license-url]  |

[license-img]: https://img.shields.io/github/license/JuliaMolSim/DftFunctionals.jl.svg?maxAge=2592000
[license-url]: https://github.com/JuliaMolSim/DftFunctionals.jl/blob/master/LICENSE

An abstract interface for evaluating DFT functionals and their derivatives
on a grid-based representation of the density.

For a small number of functionals a pure-Julia implementation is also provided,
where the focus is on readable code and a generic implementation. These functions
have not yet been tested thoroughly and likely contain numerical issues, so use them
at your own risk. Of course feel free to report or fix bugs if you encounter them!

A more stable (but less generic) implementation of DFT functionals is available
via the [Libxc.jl](https://github.com/JuliaMolSim/Libxc.jl) package.

## Packages using DftFunctionals.jl
The following packages currently use of this interface.
Please feel free to send a PR to add your package to this list.

- [DFTK.jl](https://github.com/JuliaMolSim/DFTK.jl): Plane-wave density-functional theory code in Julia
- [Libxc.jl](https://github.com/JuliaMolSim/Libxc.jl)(via an [experimental interface in DFTK](https://github.com/JuliaMolSim/DFTK.jl/tree/master/src/DispatchFunctional.jl))

```@meta
CurrentModule = DftFunctionals
```

# DftFunctionals

DftFunctionals provides an [Interface](@ref) for evaluating
exchange-correlation functionals for density-functional theory calculations.

For a small number of functionals a [Generic functional implementation](@ref)
in pure Julia is also provided,
where the focus for the moment is on readable code rather than speed.

The package is at an early stage, so please [report bugs](https://github.com/JuliaMolSim/DftFunctionals.jl/issues) or numerical issues if you find them.

A more stable (but less generic) implementation of DFT functionals is available
via the [Libxc.jl](https://github.com/JuliaMolSim/Libxc.jl) package.
Currently Libxc is not yet compatible with DftFunctionals out of the box
(but there exists an [experimental wrapper](https://github.com/JuliaMolSim/DFTK.jl/tree/master/src/DispatchFunctional.jl) in [DFTK.jl](https://dftk.org)).
This will change in the future.

## Usage
Input data to DftFunctionals is the density ``ρ`` (and its derivatives such as ``σ = \nabla \rho \cdot \nabla \rho`` etc.) on a real-space grid.
This data needs to be passed as a `(n_spin, n_p)` matrix,
where `n_spin` is the number of spin components and `n_p` is the number of
data points. Notice that three-dimensional density therefore need to be flattened
before calling anything from DftFunctionals.jl. For more details on the shapes
of the input and output arrays, see the documentation of [`potential_terms`](@ref).

```@example
using DftFunctionals

# Setup dummy input data
ρ = reshape([0.1, 0.2, 0.3, 0.4, 0.5], 1, :)
σ = reshape([0.2, 0.3, 0.4, 0.5, 0.6], 1, :)

# Setup the functional to use
xc = DftFunctional(:gga_x_pbe)

# Get the potential terms
terms = potential_terms(xc, ρ, σ)

# Energy per unit volume on each grid point
@show terms.e

# Potential Vσ = de / dρ on each grid point
@show terms.Vρ

nothing
```

"""
A generic DFT functional implementation. Valid `identifiers` are the
ones supported by [Libxc](https://tddft.org/programs/libxc) as symbols,
e.g. `:lda_x`, `hyb_gga_xc_b3lyp`. Only a subset of functionals is currently
available.
"""
function DftFunctional(identifier::Symbol, args...; kwargs...)
    DftFunctional(Val(identifier), args...; kwargs...)
end

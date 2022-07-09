struct DftFunctional{Identifier,Family,Kind} <: Functional{Family,Kind}
end
identifier(::DftFunctional{I}) where {I} = I


function DftFunctional{Identifier}() where {Identifier}
    id = string(Identifier)
    # TODO What about hybrids -> should not cause a different family,
    #      but it needs to be flagged somehow that this is a hybrid, e.g. as a field.
    @assert !startswith(id, "hyb_")  # Hybrids unsupported in fallbacks
    family, kind = split(id, '_')

    Family = Symbol(family)
    Kind   = Symbol(kind)
    @assert Family in (:lda, :gga)   # :mgga and :mggal unsupported in fallbacks
    @assert Kind   in (:x, :c, :xc)  # :k currently unsupported in fallbacks
    DftFunctional{Identifier,Family,Kind}()
end

"""
A generic DFT functional implementation. Valid identifiers `Id` are the
ones supported by [Libxc](https://tddft.org/programs/libxc) as symbols,
e.g. `:lda_x`, `hyb_gga_xc_b3lyp`.
"""
DftFunctional(Id::Symbol, args...; kwargs...) = DftFunctional{Id}(args...; kwargs...)

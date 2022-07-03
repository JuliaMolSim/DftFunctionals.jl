struct DftFunctional{Identifier,Family,Kind} <: Functional{Family,Kind}
end
identifier(::DftFunctional{I}) where {I} = I


"""A generic DFT functional implementation."""
function DftFunctional{Identifier}(args...; kwargs...) where {Identifier}
    family, kind = split(string(Identifier), '_')
    Family = Symbol(family)
    Kind   = Symbol(kind)
    @assert Family in (:lda, :gga, :mgga)
    @assert Kind   in (:x, :c, :xc, :k)
    DftFunctional{Identifier,Family,Kind}()
end
DftFunctional(Id::Symbol, args...; kwargs...) = DftFunctional{Id}(args...; kwargs...)

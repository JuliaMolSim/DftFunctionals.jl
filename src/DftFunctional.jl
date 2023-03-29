"""
    DftFunctional(identifier::Symbol, args...; kwargs...)

A generic DFT functional implementation. Valid `identifiers` are the
ones supported by [Libxc](https://tddft.org/programs/libxc) as symbols,
e.g. `:lda_x`, `hyb_gga_xc_b3lyp`. Only a subset of functionals is currently
available. Additional arguments and kwargs are passed to the functional
(e.g. to modify functional parameters).
"""
function DftFunctional(identifier::Symbol, args...; kwargs...)
    DftFunctional(Val(identifier), args...; kwargs...)
end

"""
    DftFunctional(T, identifier::Symbol)

Construct a functional with default setup and parameters and enforcing
the parameters to be represented using the number type T. Since parameters
are usually represented in Float64 precision, this can be used to enable
functional evaluations in Float32 or lower precision.
"""
function DftFunctional(T::Type{<:Number}, identifier::Symbol)
    fun = DftFunctional(Val(identifier))
    if T == parameter_type(fun)
        fun
    else
        change_parameters(fun, T.(parameters(fun)); keep_identifier=true)
    end
end

using GPUArraysCore: AbstractGPUArray

# There is no GPU implementation for DftFunctionals, instead we transfer
# the input GPU arrays to the CPU, evaluate the functional on the CPU, and
# copy the results back to the GPU (same architecture as input ρ).
to_cpu(x::AbstractGPUArray) = Array(x)
to_cpu(x) = x

# In-place transfer of CPU arrays `x` to an existing GPU array `y`.
# This does not allocate a new array.
to_gpu!(y::AbstractGPUArray, x::AbstractArray) = (copyto!(y, x); y)
to_gpu!(y, x) = y  # non-arrays: nothing to transfer
to_gpu!(nt_y::NamedTuple, nt_x::NamedTuple) = map(to_gpu!, nt_y, nt_x)

function cpu_fallback(cpu_fun, func, ρ::AbstractGPUArray, args...)
    args_cpu = map(to_cpu, args)
    result_cpu = cpu_fun(func, to_cpu(ρ), args_cpu...)
    result_gpu = map(result_cpu) do x
        x isa AbstractArray ? similar(ρ, eltype(x), size(x)) : x
    end
    to_gpu!(result_gpu, result_cpu)
end

# LDA
for fun in (:energy_density, :potential_terms, :kernel_terms)
    @eval function $fun(func::Functional{:lda}, ρ::AbstractGPUArray)
        cpu_fallback($fun, func, ρ)
    end
end

# GGA
for fun in (:energy_density, :potential_terms, :kernel_terms)
    @eval function $fun(func::Functional{:gga}, ρ::AbstractGPUArray, σ::AbstractGPUArray)
        cpu_fallback($fun, func, ρ, σ)
    end
end

# TODO: Add mGGA fallbacks here once mGGA functionals are implemented.

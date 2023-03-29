@doc raw"""
In the functional interface we explicitly use the symmetry ``σ_{αβ} = σ_{βα}`` for
the contracted density gradient ``σ_ij = ∇ρ_i ⋅ ∇ρ_j`` where ``i,j ∈ \{α, β\}``.
Instead of treating ``σ`` to have the spin components ``σ_{αα}``, ``σ_{αβ}``, ``σ_{βα}``
and ``σ_{ββ}``, we consider it to have three pseudo-spin components ``σ_{αα}``, ``σ_x``
and ``σ_{ββ}``, where ``σ_x = (σ_{αβ} + σ_{βα})/2 = σ_{αβ}``. Input arrays like `σ` or
output arrays like `Vσ` or `Vρσ` will therefore feature a spin axis of length 3 refering
to the ``σ_{αα}``, ``σ_{αβ}`` and ``σ_{ββ}`` components / derivatives respectively.
E.g. `σ` is of shape `(3, n_p)` (where `n_p` is the number of points where evaluation takes
place.

This function maps the "unfolded" spin tuple `(i, j)` in ``σ_ij`` to the corresponding
place in the length-3 spin axis of `σ`.
"""
function spinindex_σ(s::Int, t::Int)
    s == 1 && t == 1 && return 1
    s == 2 && t == 2 && return 3
    return 2
end


"""
Determine the working precision of the functional evaluation. This takes the type
of the parameters of Functional into account, but only in case this type is not a plain
number (e.g. to allow duals to get pushed through).
"""
function working_precision(func::Functional, T, S...)
    working_precision_(eltype(parameters(func)), T, S...)
end
working_precision_(::Type{<:Number}, T::Type, S...) = promote_type(T, S...)
working_precision_(PT::Type, T, S...)               = promote_type(PT, T, S...)

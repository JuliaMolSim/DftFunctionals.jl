using DftFunctionals
using ComponentArrays
using Test
include("libxc.jl")

@testset "DftFunctionals.jl" begin
    gga_fallback = (
        :gga_x_pbe,     :gga_c_pbe,     :gga_x_pbe_r,
        :gga_x_pbe_sol, :gga_c_pbe_sol,
        :gga_x_pbefe,   :gga_c_pbefe,
        :gga_x_xpbe,    :gga_c_xpbe,
        :gga_x_pbe_mol, :gga_c_pbe_mol,
        :gga_x_apbe,    :gga_c_apbe,
    )

@testset "LDA Functional construction" begin
    let f = DftFunctional(:lda_x)
        @test kind(f)       == :x
        @test family(f)     == :lda
        @test identifier(f) == :lda_x
        @test !needs_σ(f)
        @test !needs_τ(f)
        @test !needs_Δρ(f)
        @test has_energy(f)
    end

    for id in (:lda_c_vwn, :lda_c_pw)
        f = DftFunctional(id)
        @test kind(f)       == :c
        @test family(f)     == :lda
        @test identifier(f) == id
    end
end

@testset "GGA Functional construction" begin
    for id in gga_fallback
        f = DftFunctional(id)
        @test family(f) == :gga
        @test identifier(f) == id
        @test needs_σ(f)
        @test !needs_τ(f)
        @test !needs_Δρ(f)
    end
end


@testset "LDA potential (without spin)" begin
    ρ     = [0.1, 0.2, 0.3, 0.4, 0.5]
    ref_e = [-0.0342808612, -0.0863823573, -0.1483246721, -0.2176699007, -0.2930972406]
    ref_v = [-0.4570781497, -0.5758823823, -0.6592207650, -0.7255663357, -0.7815926418]

    ρ      = reshape(ρ, 1, :)
    ref_v  = reshape(ref_v, 1, :)
    result = potential_terms(DftFunctional(:lda_x), ρ)
    @test result.e  ≈ ref_e atol=1e-7
    @test result.Vρ ≈ ref_v atol=1e-7
end

@testset "GGA potential (without spin)" begin
    ρ      = [0.1, 0.2, 0.3, 0.4, 0.5]
    σ      = [0.2, 0.3, 0.4, 0.5, 0.6]
    ref_e  = [-0.0452597523, -0.0957755740, -0.1562022650, -0.2245711188, -0.2993306248]
    ref_Vρ = [-0.4273142446, -0.5301998810, -0.6288348093, -0.7043765616, -0.7658496683]
    ref_Vσ = [-0.0330279599, -0.0270759865, -0.0183930411, -0.0132581503, -0.0101141659]

    ρ      = reshape(ρ, 1, :)
    σ      = reshape(σ, 1, :)
    ref_Vρ = reshape(ref_Vρ, 1, :)
    ref_Vσ = reshape(ref_Vσ, 1, :)
    result = potential_terms(DftFunctional(:gga_x_pbe), ρ, σ)
    @test result.e  ≈ ref_e  atol=1e-7
    @test result.Vρ ≈ ref_Vρ atol=1e-7
    @test result.Vσ ≈ ref_Vσ atol=1e-7
end


@testset "Fallback <-> Libxc (LDA, without spin)" begin
    # Note: Libxc defaults to the PW correlation functional as published,
    #       we to an improved form with better constants ... that's why we need
    #       the improved=false below.
    for func in (DftFunctional(:lda_x),
                 DftFunctional(:lda_c_vwn),
                 DftFunctional(:lda_c_pw; improved=false))
        @testset "$(identifier(func))" begin
            n_p   = 100
            ρ     = abs.(randn(1, n_p))
            εref  = similar(ρ, n_p)
            Vref  = similar(ρ)
            V2ref = similar(ρ)

            ptr = xc_functional_alloc(identifier(func))
            xc_lda(ptr, n_p, ρ, εref, Vref, V2ref, C_NULL, C_NULL)
            xc_functional_free(ptr)

            eref  = εref .* ρ[1, :]
            V2ref = reshape(V2ref, 1, 1, :)

            # Compute in fallback implementation in elevated precision
            result = kernel_terms(func, Array{BigFloat}(ρ))
            @test result.e   ≈ eref  atol=5e-13
            @test result.Vρ  ≈ Vref  atol=5e-13
            @test result.Vρρ ≈ V2ref atol=5e-13

            # Check floating-point type consistency:
            e = DftFunctionals.energy(func, rand(Float32))
            @test eltype(e) == Float32
        end
    end
end

@testset "Fallback <-> Libxc (GGA, without spin)" begin
    for func_name in gga_fallback
        @testset "$func_name" begin
            func  = DftFunctional(func_name)

            n_p    = 100
            ρ      = abs.(randn(1, n_p))
            σ      = abs.(randn(1, n_p))
            εref   = similar(ρ, n_p)
            Vρref  = similar(ρ)
            Vσref  = similar(ρ)
            Vρρref = similar(ρ)
            Vρσref = similar(ρ)
            Vσσref = similar(ρ)


            ptr = xc_functional_alloc(identifier(func))
            xc_gga(ptr, n_p, ρ, σ, εref, Vρref, Vσref, Vρρref, Vρσref, Vσσref, C_NULL,
                   C_NULL, C_NULL, C_NULL, C_NULL, C_NULL, C_NULL, C_NULL, C_NULL)
            xc_functional_free(ptr)

            eref  = εref .* ρ[1, :]
            Vρρref = reshape(Vρρref, 1, 1, :)
            Vρσref = reshape(Vρσref, 1, 1, :)
            Vσσref = reshape(Vσσref, 1, 1, :)

            # Compute in fallback implementation in elevated precision
            result = kernel_terms(func, Array{BigFloat}(ρ), Array{BigFloat}(σ))
            @test result.e    ≈ eref   atol=5e-13
            @test result.Vρ   ≈ Vρref  atol=5e-13
            @test result.Vσ   ≈ Vσref  atol=5e-13
            @test result.Vρρ  ≈ Vρρref atol=5e-13
            @test result.Vρσ  ≈ Vρσref atol=5e-13
            @test result.Vσσ  ≈ Vσσref atol=5e-13

            # Ensure Float32 evaluation works (if parameters are forcibly down-casted to 32bit)
            func32 = DftFunctional(Float32, func_name)
            e = DftFunctionals.energy(func32, rand(Float32), rand(Float32))
            @test eltype(e) == Float32
        end
    end
end

@testset "Test spinindex_σ" begin
    @test spinindex_σ(1, 1) == 1
    @test spinindex_σ(1, 2) == 2
    @test spinindex_σ(2, 1) == 2
    @test spinindex_σ(2, 2) == 3
end

@testset "PBE functionals" begin
    pbe = DftFunctional(:gga_x_pbe)
    @test :μ in keys(parameters(pbe))
    @test :κ in keys(parameters(pbe))

    pbemod = change_parameters(pbe, ComponentArray(;μ=12, κ=1.2))
    @test parameters(pbemod).μ == 12
    @test parameters(pbemod).κ == 1.2

    pbemod = change_parameters(DftFunctional(:gga_c_pbe), ComponentArray(;β=12, γ=1.2))
    @test parameters(pbemod).β == 12
    @test parameters(pbemod).γ == 1.2

    μ = rand()
    @test μ ≈ DftFunctionals.pbe_μ_from_β(DftFunctionals.pbe_β_from_μ(μ))
end

@testset "PBE parameter derivatives" begin
    using ForwardDiff

    pbe = DftFunctional(:gga_x_pbe)

    ρ = [0.1, 0.2, 0.3, 0.4, 0.5]
    σ = [0.2, 0.3, 0.4, 0.5, 0.6]
    ρ = reshape(ρ, 1, :)
    σ = reshape(σ, 1, :)

    θ = ComponentArray(; parameters(pbe)...)
    egrad = ForwardDiff.jacobian(θ) do θ
        potential_terms(change_parameters(pbe, θ), ρ, σ).e
    end

    egrad_fd = let ε=1e-5
        δ = zero(θ)
        δ[2] = ε

        (  potential_terms(change_parameters(pbe, θ + δ), ρ, σ).e
         - potential_terms(change_parameters(pbe, θ - δ), ρ, σ).e) / 2ε
    end

    @test maximum(abs, egrad[:, 2] - egrad_fd) < 1e-5
end
end

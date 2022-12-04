using Test, SpecialFunctions

include("cref/fet_cref.jl")
include("../src/fet.jl")

@testset "Fisher's Exact Test" begin
    @testset "kf_lgamma" begin
        kf_lgamma = fet.kf_lgamma
        @test kf_lgamma_cref(1.0) ≈ kf_lgamma(1.0)
        @test kf_lgamma_cref(3.1415926) ≈ kf_lgamma(3.1415926)
        for i in 1:9
            n = 10^i + rand()
            @test kf_lgamma_cref(n) ≈ kf_lgamma(n)
        end
    end

    @testset "kf_erfc" begin
        kf_erfc = fet.kf_erfc
        @test kf_erfc_cref(0.0) ≈ kf_erfc(0.0)
        @test kf_erfc_cref(2.0) ≈ kf_erfc(2.0)

        for i in 1:9
            n = 10^i + rand()
            @test kf_erfc_cref(n) ≈ kf_erfc(n)
        end
    end

    @testset "_kf_gammap" begin
        _kf_gammap = fet._kf_gammap
        @test _kf_gammap_cref(2.0, 3.0) ≈ _kf_gammap(2.0, 3.0)
        @test _kf_gammap_cref(0.0, 0.1) ≈ _kf_gammap(0.0, 0.1)
        @test _kf_gammap_cref(1.0, 0.1) ≈ _kf_gammap(1.0, 0.1)
    end

    @testset "_kf_gammaq" begin
        _kf_gammaq = fet._kf_gammaq
        @test _kf_gammaq_cref(2.0, 3.0) ≈ _kf_gammaq(2.0, 3.0)
        @test _kf_gammaq_cref(0.0, 0.1) ≈ _kf_gammaq(0.0, 0.1)
        @test _kf_gammaq_cref(1.0, 0.1) ≈ _kf_gammaq(1.0, 0.1)
    end

    @testset "kf_gammap" begin
        kf_gammap = fet.kf_gammap
        @test kf_gammap_cref(2.0, 3.0) ≈ kf_gammap(2.0, 3.0)
        @test kf_gammap_cref(0.0, 0.1) ≈ kf_gammap(0.0, 0.1)
        @test kf_gammap_cref(1.0, 0.1) ≈ kf_gammap(1.0, 0.1)
    end

    @testset "kf_gammaq" begin
        kf_gammaq = fet.kf_gammaq
        @test kf_gammaq_cref(2.0, 3.0) ≈ kf_gammaq(2.0, 3.0)
        @test kf_gammaq_cref(0.0, 0.1) ≈ kf_gammaq(0.0, 0.1)
        @test kf_gammaq_cref(1.0, 0.1) ≈ kf_gammaq(1.0, 0.1)
    end

    @testset "kf_betai_aux" begin
        kf_betai_aux = fet.kf_betai_aux
        a, b, x = 0.1, 0.1, 0.1
        @test kf_betai_aux_cref(a, b, x) ≈ kf_betai_aux(a, b, x)
    end

    # @testset "kf_betai" begin
        # kf_betai = fet.kf_betai
        # a, b, x = 0.1, 0.8, 0.4
        # @test kf_betai_cref(a, b, x) ≈ kf_betai(a, b, x)
    # end

    @testset "lbinom" begin
        lbinom = fet.lbinom
        n, k = 0, 114
        @test lbinom_cref(n, k) ≈ lbinom(n, k)
    end

    @testset "hypergeo" begin
        hypergeo = fet.hypergeo
        n11, n1_, n_1, n = 1656, 32106, 17466, 251573
        @test hypergeo_cref(n11, n1_, n_1, n) ≈ hypergeo(n11, n1_, n_1, n)

        n11, n1_, n_1, n = 3, 4, 4, 8
        @test hypergeo_cref(n11, n1_, n_1, n) ≈ hypergeo(n11, n1_, n_1, n)
    end

    @testset "hypergeo_acc" begin
        hypergeo_acc = fet.hypergeo_acc

        n11, n1_, n_1, n = 3, 0, 0, 0
        aux = hgacc_t(4, 4, 4, 8, 0.014286)
        aux_ref = Ptr{hgacc_t}(pointer_from_objref(aux))

        ref = hypergeo_acc_cref(n11, n1_, n_1, n, aux_ref)

        aux = hgacc_t(4, 4, 4, 8, 0.014286)
        aux_ref = Ptr{hgacc_t}(pointer_from_objref(aux))
        val = hypergeo_acc(n11, n1_, n_1, n, aux)
        @test ref ≈ val
    end

    @testset "kt_fisher_exact" begin
        kt_fisher_exact = fet.kt_fisher_exact

        # @test kt_fisher_exact_cref(1656, 32106, 17466, 251573)[4] ≈ 1.825331208626037*10^(-31)
        # @test kt_fisher_exact_cref(3, 1, 1, 3)[4] ≈ 0.48571429
        # @test kt_fisher_exact(1656, 32106, 17466, 251573)[4] ≈ 1.825331208626037*10^(-31)
        # @test kt_fisher_exact(3, 1, 1, 3)[4] ≈ 0.48571429
    end
end;

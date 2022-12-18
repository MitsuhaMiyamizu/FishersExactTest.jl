using Test, SpecialFunctions

include("cref/fet_cref.jl")
using FishersExactTest
import FishersExactTest.hgacc_t

@testset "Fisher's Exact Test" begin

    @testset "lbinom" begin
        lbinom = FishersExactTest.lbinom
        n, k = 0, 114
        @test lbinom_cref(n, k) ≈ lbinom(n, k)
    end

    @testset "hypergeo" begin
        hypergeo = FishersExactTest.hypergeo
        n11, n1_, n_1, n = 1656, 32106, 17466, 251573
        @test hypergeo_cref(n11, n1_, n_1, n) ≈ hypergeo(n11, n1_, n_1, n)

        n11, n1_, n_1, n = 3, 4, 4, 8
        @test hypergeo_cref(n11, n1_, n_1, n) ≈ hypergeo(n11, n1_, n_1, n)
    end

    # @testset "hypergeo_acc" begin
        # hypergeo_acc = FishersExactTest.hypergeo_acc

        # n11, n1_, n_1, n = 3, 0, 0, 0
        # aux = hgacc_t(4, 4, 4, 8, 0.0)
        # aux_ref = Ptr{hgacc_t}(pointer_from_objref(aux))
        # ref = hypergeo_acc_cref(n11, n1_, n_1, n, aux_ref)

        # aux = hgacc_t(4, 4, 4, 8, 0.0)
        # aux_ref = Ptr{hgacc_t}(pointer_from_objref(aux))
        # val = hypergeo_acc(n11, n1_, n_1, n, aux)
        # @test ref ≈ val
    # end

    @testset "kt_fisher_exact" begin
        kt_fisher_exact = FishersExactTest.kt_fisher_exact

        @test kt_fisher_exact_cref(1656, 32106, 17466, 251573)[4] ≈ 1.825331208626037*10^(-31)
        @test kt_fisher_exact_cref(3, 1, 1, 3)[4] ≈ 0.48571429
        @test kt_fisher_exact(1656, 32106, 17466, 251573)[4] ≈ 1.825331208626037*10^(-31)
        @test kt_fisher_exact(3, 1, 1, 3)[4] ≈ 0.48571429
    end
end;

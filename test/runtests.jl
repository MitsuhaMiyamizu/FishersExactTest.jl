using fet
using Test

include("cref/fet_cref.jl")

@testset "test suit for Fisher's exact test" begin
           @test fisher_exact_cref(1656, 32106, 17466, 251573)[4] ≈ 1.825331208626037*10^(-31)
           @test isapprox(fisher_exact_cref(3, 1, 1, 3)[4], 0.48571429)
end;

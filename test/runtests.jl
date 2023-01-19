using Colocalization
using Test
using Random

@testset "Colocalization.jl" begin
    @testset "Coloc" begin
        X, Y = 100, 200
        Random.seed!(42)
        for i in 1:100
            xs = Random.rand(X, Y)
            ys = Random.rand(X, Y) .+ .2
            res = colocalize_all(xs, ys)
            for k in keys(Colocalization.metrics)
                μ, md, σ, min, max, q1, q3, q95, q99, nans = describe_array(res[k])
                @test μ >= 0
                @test md >= 0
                @test σ >= 0
                @test min <= max
                @test q1 <= q3
                @test q95 <= q99
                @test 0 <= nans < X*Y
            end
            df = summarize_colocalization(res, "f1", "f2")
            @test size(df, 1) == length(Colocalization.metrics)
            @test !isnothing(df)
        end
    end
end

using Colocalization
using Test
using Random
using Logging

@testset "Colocalization.jl" begin
    @testset "Coloc" begin
        c = global_logger()
        global_logger(NullLogger())
        X, Y = 100, 200
        Random.seed!(42)
        for i in 1:100
            xs = Random.rand(X, Y)
            ys = Random.rand(X, Y) .+ .2
            res = colocalize_all(xs, ys)
            for k in keys(Colocalization.metrics)
                μ, md, σ, m, M, q1, q3, q95, q99, nans = describe_array(res[k])
                @test μ >= -1
                @test md >= -1
                @test σ >= 0
                @test m <= M
                @test q1 <= q3
                @test q95 <= q99
                @test 0 <= nans < X*Y
            end
            df = summarize_colocalization(res, "f1", "f2")
            @test size(df, 1) == length(Colocalization.metrics)
            @test !isnothing(df)
        end
        global_logger(NullLogger())
    end
    @testset "tomask" begin
        Random.seed!(42)
        X = 10
        Y = 20
        Z = 2
        img1 = rand(X,Y,Z)
        i2 = copy(img1)
        msk = tomask(img1)
        @test all(i2 .== img1)
        @test all(msk[i2 .>0].== 1)
        @test all(msk[i2 .== 0].== 0)
    end

    @testset "nowindow" begin
        X = ones(20, 20)
        Y = ones(20, 20)
        r=colocalize_nowindow(X, Y)
        @test iszero(r)
        X = ones(20, 20) .+ rand(20, 20)./100
        Y = ones(20, 20) .- rand(20, 20)./100
        r=colocalize_nowindow(X, Y)
        @test 0 < abs(r) <= .2
    end

    @testset "aszero" begin
        Random.seed!(42)
        X = 10
        Y = 20
        Z = 2
        img1 = rand(X,Y,Z)
        img2 = copy(img1)
        img0 = aszero(img1)
        @test all(img0 .== zero(eltype(img1)))
        @test all(img1 .== img2)
    end
    @testset "full_coloc" begin
        X = zeros(10, 10)
        Y = zeros(10, 10)
        UM = union_mask(X, Y)
        @test sum(UM) == 0
        @test sum(intersection_mask(X, Y)) == 0
        X[4:6, 4:6].=1
        Y[5:7, 5:7].=1
        @test sum(union_mask(X, Y)) == 14
        @test sum(intersection_mask(X, Y)) == 4
    end

    @testset "df" begin
        X = zeros(20, 20)
        Y = zeros(20, 20)
        X[4:6, 4:6].=1
        Y[5:7, 5:7].=1
        X[10:12, 10:12] .= 1
        Y[15:16, 15:16] .= 1
        xr, yr = union_distance_mask(X, Y, 8)
        @test sum(xr) == sum(X)
        @test sum(yr) == sum(Y)
    end
end

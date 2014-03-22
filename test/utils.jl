using Matfacgrf
using Base.Test

import Matfacgrf.residual

function testResidual()
    m, k, n = 10, 5, 12
    A = sprand(m, n, .3)
    W = rand(m, k)
    H = rand(k, n)
    cr = residual(A, W, H, 2)
    rr = residual(A, W, H, 1)
    tr = residual(A, W, H)
    @test_approx_eq sum(cr) tr# 10e-10
    @test_approx_eq sum(rr) tr# 10e-10
    @test_approx_eq sum(rr) sum(cr)
    @show cr, rr, tr
    info("summing both dimensions gives same answer")

    B = sparse(W * H)
    @show Matfacgrf.residual(B,W,H)
    @test_approx_eq sum(Matfacgrf.residual(B, W, H)) 0.0    #10e-10
    @test_approx_eq sum(Matfacgrf.residual(B, W, H, 2)) 0.0 #10e-10
    @test_approx_eq sum(Matfacgrf.residual(B, W, H, 1)) 0.0 #10e-10
    #test throws
    @test_throws Matfacgrf.residual(A, W, [1 2])
    info("test of residual computations passed.")
end

function testDenseResidual()
    m, k, n = 10, 5, 12
    W = rand(m,k)
    H = rand(k,n)
    A = W*H
    tr = residual(A, W, H)
    @test tr == 0.0
    info("Residual works of dense matrices too")
end

testResidual()
testDenseResidual()

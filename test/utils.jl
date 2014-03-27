using Matfacgrf
using Base.Test

import Matfacgrf.residual
tolerance = 10e-10
function testResidual()
    m, k, n = 10, 5, 12
    A = sprand(m, n, .3)
    W = rand(m, k)
    H = rand(k, n)
    cr = residual(A, W, H, 2)
    rr = residual(A, W, H, 1)
    tr = residual(A, W, H)
    @test_approx_eq_eps sum(cr) tr tolerance
    @test_approx_eq_eps sum(rr) tr tolerance
    @test_approx_eq_eps sum(rr) sum(cr) tolerance
    @show cr, rr, tr
    info("summing both dimensions gives same answer")

    B = sparse(W * H)
    @show Matfacgrf.residual(B,W,H)
    @test_approx_eq_eps sum(Matfacgrf.residual(B, W, H)) 0.0    tolerance
    @test_approx_eq_eps sum(Matfacgrf.residual(B, W, H, 2)) 0.0 tolerance
    @test_approx_eq_eps sum(Matfacgrf.residual(B, W, H, 1)) 0.0 tolerance
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

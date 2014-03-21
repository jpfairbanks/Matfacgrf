using Base.Test
function testResidual()
    m, k, n = 10, 5, 12
    A = sprand(m,n, .3)
    W = rand(m,k)
    H = rand(k,n)
    cr = residual(A, W, H, 2)
    rr = residual(A, W, H, 1)
    tr = residual(A, W, H)
    @test sum(cr) == tr
    @test sum(rr) == tr
    @test sum(rr) == sum(cr)

    B = sparse(W*H)
    @test residual(B,W,H) == 0.0
    #test throws
    @test_throws residual(A, W, [1 2])
end

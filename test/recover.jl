using Base.Test
using NMF
using Matfacgrf

test_tolerance = 10e-12
function denseWH(m,k,n)
    #make a W,H pair that is random for testing purposes.
    W = rand(m, k)
    H = rand(k, n)
    return W,H
end

function eyeWH(m,k,n)
    #make a W,H pair that is simple for testing purposes.
    W = eye(m,k)
    H = eye(k,n)
    return W,H
end

function sparseWH(m,k,n,densityW, densityH)
    #make a random sparse W,H pair for testing.
    W = full(sprand(m, k, densityW))
    H = full(sprand(k, n, densityH))
    return W,H
end

function testInitialization()
    info("testInitialization of Hals algorithm and matrices")
    m = 5
    n = 4
    k = 3
    alg = Matfacgrf.HierarchicalALS(
        maxiter = 50,
        tol = 10e-6,
        lambda = 0,
        verbose = false
        )
    W,H = Matfacgrf.randinit(m,n,k;normalize=true, zeroh=true)
    info("test bad initialization produces error")
    @test_throws Matfacgrf.solve!(alg, rand(m,n), W, H)
    W,H = Matfacgrf.randinit(m,n,k;normalize=true, zeroh=false)
    info("test return type")
    @test typeof(Matfacgrf.solve!(alg, rand(m,n), W, H)) <: NMF.Result
    W = zeros(m,k)
    info("test sparse call")
    @test Matfacgrf.solve!(alg, sprand(m, n, 0.3), W, H) <: NMF.Result
end

function finderrors(W,H,k)
    A = W*H
    tol = 10e-14
    numIter = 1000
    beta = 0.0
    factors = Matfacgrf.hals(A,W,H,k,tol,numIter,beta, false)
    newW = factors.W
    newH = factors.H
    werr = vecnorm(W - newW)/vecnorm(W)
    herr = vecnorm(H - newH)/vecnorm(H)
    residual = A - (newW * newH)
    residualerr = vecnorm(residual)/vecnorm(A)
    return werr, herr, residualerr
end
function testRecover()
    info("testRecover: tests that we can recover the factors for simple cases.")
    #the goal is to draw two random matrices H and W and make sure that
    #we can recover them with hals.m
    n = 60
    m = 70
    k = 50
    info("dense test")
    W, H = denseWH(m,k,n)
    @test_approx_eq_eps finderrors(W,H,k)[3] 0.0 test_tolerance
    info("sparse test")
    W, H = sparseWH(m,k,n,.3,.4)
    @test_approx_eq_eps finderrors(W,H,k)[3] 0.0 test_tolerance
    info("eye test")
    W, H = eyeWH(k+1,k,k)
    #check that H==0 throws an error
    @test_throws finderrors(W,H,k)
    W, H = eyeWH(k,k,k)
    @test_approx_eq_eps finderrors(W,H,k)[3] 0.0 test_tolerance
end

testInitialization()
testRecover()

using Base.Test
using NMF

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
    m = 5
    n = 4
    k = 3
    tol = 10e-6
    numIter = 50
    beta = 0
    W,H = randinit(m,n,k;normalize=true, zeroh=true)
    info("first test")
    @test_throws hals(rand(m,n), W, H, k, tol, numIter, beta)
    W,H = randinit(m,n,k;normalize=true, zeroh=false)
    info("Second test")
    @test typeof(hals(rand(m,n), W, H, k, tol, numIter, beta)) <: NMF.Result
    W = zeros(m,k)
    info("Third Test")
    @show hals(rand(m,n), W, H, k, tol, numIter, beta)
end

function finderrors(W,H,k)
    A = W*H
    tol = 10e-14
    numIter = 1000
    beta = 0.0
    factors = hals(A,W,H,k,tol,numIter,beta)
    newW = factors.W
    newH = factors.H
    werr = normfro(W - newW)/normfro(W)
    herr = normfro(H - newH)/normfro(H)
    residual = A - (newW * newH)
    residualerr = normfro(residual)/normfro(A)
    return werr, herr, residualerr
end
function testRecover()
    #the goal is to draw two random matrices H and W and make sure that
    #we can recover them with hals.m
    n = 600
    m = 700
    k = 500
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

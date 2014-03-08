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
function finderrors(W,H,k)
    A = W*H
    tol = .00001
    numIter = 100
    beta = 0.0001
    (newW, newH, errlog) = hals(A,W,H,k,tol,numIter,beta)
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
    @show Case = "dense test"
    W, H = denseWH(m,k,n)
    @show finderrors(W,H,k)
    @show Case = "sparse test"
    W, H = sparseWH(m,k,n,.3,.4)
    @show finderrors(W,H,k)
    @show Case = "eye test"
    W, H = eyeWH(m,k,n)
    @show finderrors(W,H,k)
end

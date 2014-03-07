%function testRecover()
    %the goal is to draw two random matrices H and W and make sure that
    %we can recover them with hals.m
    n = 6;
    m = 7;
    k = 5;
    %eps=0.01
    density = .3;
    W = sprand(m,k,density)
    H = sprand(k,n,density)
    A = W*H
    tol = .0001;
    numIter = 50;
    beta = 1;
    [newW, newH, errlog] = hals(A,W,H,k,tol,numIter,beta)
    werr = norm(W - newW,'fro')/norm(W, 'fro')
    herr = norm(H - newH,'fro')/norm(H, 'fro')
    residual = A - (newW * newH)
%end
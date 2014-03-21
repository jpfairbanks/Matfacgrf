function [werr, herr, residualerr] = finderrors(W,H)
    A = W*H;
    tol = .00001;
    numIter = 100;
    beta = 0.0001;
    [newW, newH, errlog] = hals(A,W,H,k,tol,numIter,beta);
    werr = norm(W - newW,'fro')/norm(W, 'fro');
    herr = norm(H - newH,'fro')/norm(H, 'fro');
    residual = A - (newW * newH);
    residualerr = norm(residual, 'fro')/norm(A,'fro');
end

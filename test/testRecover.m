function testRecover()
    %the goal is to draw two random matrices H and W and make sure that
    %we can recover them with hals.m
    n = 6;
    m = 7;
    k = 5;
    Case = 'dense test'
    W, H = denseWH(m,k,n)
    A = W * H
    finderrors(W,H)
    Case = 'sparse test'
    W, H = sparseWH(m,k,n)
    A = W * H
    finderrors(W,H)
    Case = 'eye test'
    W, H = eyeWH(m,k,n)
    A = W * H
    finderrors(W,H)
end

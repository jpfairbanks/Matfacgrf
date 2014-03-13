function hals(A,Winit,Hinit,k,tolerance,numIterations,beta=0)
#hals stands for Hierarchical Alternating Least Squares
#solves A=WH
#A = mxn matrix
#W = mxk
#H = kxn
#k = rank of approximation.
#implementation of the algorithm2 from
#http://www.bsp.brain.riken.jp/publications/2009/Cichocki-Phan-IEICE_col.pdf
#All entries of A must be nonnegative.
#A must not contain any rows with no zero elements. We leave the user to decide how to ensure this.
assert all(sum(A,1).>0) && all(sum(A,2).>0)
W=Winit
H=Hinit
epsilon = eps(Float64)
prevError=normfro(A-W*H)
currError = prevError+1
currentIteration=1
errChange=zeros(1,numIterations)
while (abs(currError-prevError)>tolerance && currentIteration<numIterations)
    #update W
    AHt=A*H'
    HHt=H*H'
    #to avoid divide by zero error.
    HHtDiag=diag(HHt)
    @assert all(HHtDiag.>=0)
    for x=1:k
        Wx = W[:,x] + (AHt[:,x]-W*HHt[:,x])/HHtDiag[x]
        Wx[Wx.<epsilon]=epsilon
        W[:,x]=Wx
    end
    #update H
    WtA=W'*A
    WtW=W'*W
    #to avoid divide by zero error.
    WtWDiag=diag(WtW)
    @assert all(WtWDiag.>=0)
    #temporary variable prevH to overcome parfor error.
    for x=1:k
        Hx = H[x,:]+(WtA[x,:]-WtW[x,:]*H)/WtWDiag[x]
        Hx=Hx-beta/WtWDiag[x]
        Hx[Hx.<epsilon]=epsilon
        H[x,:]=Hx
    end

    if (currentIteration>1)
        prevError=currError
    end
    errChange[currentIteration]=prevError
    currError=normfro(A-W*H)
    currentIteration=currentIteration+1
end
return  W, H, errChange
end

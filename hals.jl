function hals(A,Winit,Hinit,k,tolerance,numIterations,beta)
#hals stands for Hierarchical Alternating Least Squares
#solves A=WH
#A = mxn matrix
#W = mxk
#H = kxn
#k = rank of approximation.
#implementation of the algorithm2 from
#http://www.bsp.brain.riken.jp/publications/2009/Cichocki-Phan-IEICE_col.pdf
W=Winit
H=Hinit
epsilon = eps(Float64)
prevError=normfro(A-W*H)
currError = prevError+1
currentIteration=1
errChange=zeros(1,numIterations)
while (abs(currError-prevError)>tolerance && currentIteration<numIterations)
    #update W
    #temporary variable prevW to overcome parfor error.
    AHt=A*H'
    HHt=H*H'
    #to avoid divide by zero error.
    HHtDiag=diag(HHt)
    HHtDiag[HHtDiag.==0]=epsilon
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
    WtWDiag[WtWDiag.==0]=epsilon
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

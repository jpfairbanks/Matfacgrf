function [ W,H,errChange ] = hals( A,Winit,Hinit,k,tolerance,numIterations)
%hals stands for Hierarchical Alternating Least Squares
%both tolerance and beta must be very small in order to obtain good results
%solves A=WH
%A = mxn matrix
%W = mxk
%H = kxn
%k = rank of approximation.
%implementation of the algorithm2 from
%http://www.bsp.brain.riken.jp/publications/2009/Cichocki-Phan-IEICE_col.pdf
W=Winit;
H=Hinit; %Hinit must not result in zero diagonal for HH'
prevError=norm(A-W*H,'fro');
currError = prevError+1;
currentIteration=1;
errChange=zeros(1,numIterations);
while (abs(currError-prevError)>tolerance && currentIteration<numIterations)
    %update W;
    %temporary variable prevW to overcome parfor error.
    AHt=A*H';
    HHt=H*H';
    %assert divide by zero error.
    HHtDiag=diag(HHt);
    assert(isempty(HHtDiag==0));
    for x=1:k
        Wx = W(:,x) + (AHt(:,x)-W*HHt(:,x))/HHtDiag(x);
        Wx(Wx<eps)=eps;
        W(:,x)=Wx;
    end
    %update H
    WtA=W'*A;
    WtW=W'*W;
    %to avoid divide by zero error.
    WtWDiag=diag(WtW);
    assert(isempty(WtWDiag==0));
    for x=1:k
        Hx = H(x,:)+(WtA(x,:)-WtW(x,:)*H)/WtWDiag(x);
        Hx(Hx<eps)=eps;
        H(x,:)=Hx;
    end
    if (currentIteration>1)
        prevError=currError;
    end
    errChange(currentIteration)=prevError;
    currError=norm(A-W*H,'fro');
    currentIteration=currentIteration+1;
end
end

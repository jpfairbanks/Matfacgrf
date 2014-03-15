using NMF
function hals(A, Winit, Hinit, k, tolerance, maxiter, beta=0)
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
    @assert all(sum(A,1).>0) && all(sum(A,2).>0)
    W=Winit
    H=Hinit
    epsilon = eps(Float64)
    prevError=normfro(A-W*H)
    currError = prevError+1
    currentIteration=1
    errChange=zeros(1,maxiter)
    converged = false
    while (!converged && currentIteration<maxiter)

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
        for x=1:k
            Hx = H[x,:]+(WtA[x,:]-WtW[x,:]*H)/WtWDiag[x]
            Hx=Hx-beta/WtWDiag[x]
            Hx[Hx.<epsilon]=epsilon
            H[x,:]=Hx
        end

        #check convergence
        if (currentIteration>1)
            prevError=currError
        end
        errChange[currentIteration]=prevError
        currError=normfro(A-W*H)
        converged = abs(currError-prevError)>tolerance
        currentIteration=currentIteration+1
    end
    return NMF.Result(W, H, currentIteration, converged, currError)
end

function randinit(m::Integer,n::Integer, k::Integer; normalize::Bool=false, zeroh::Bool=false)
   #X is m by n and we want a rank k factorization.
   W = rand(m, k)
   if normalize
       normalize1_cols!(W)
   end
   H = zeroh ? zeros(k, n) : rand(k, n)
   return (W, H)::(Matrix{Float64}, Matrix{Float64})
end


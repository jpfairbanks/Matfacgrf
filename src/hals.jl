#
# Hierarchical Alternating Least Squares Algorithm for NMF.
#
# Reference:
# @article{cichocki2009fast,
#   title={Fast local algorithms for large scale nonnegative matrix and tensor factorizations},
#   author={Cichocki, Andrzej and Anh-Huy, PHAN},
#   journal={IEICE transactions on fundamentals of electronics, communications and computer sciences},
#   volume={92},
#   number={3},
#   pages={708--721},
#   year={2009},
#   publisher={The Institute of Electronics, Information and Communication Engineers},
#   url={http://www.bsp.brain.riken.jp/publications/2009/Cichocki-Phan-IEICE_col.pdf}
# }
#
# original MATLAB implementation by Ramakrishnan Kannan at Georgia Tech.
# original Julia implementation by James Fairbanks at Georgia Tech.

using NMF
type HierarchicalALS
    maxiter::Int    # maximum number of iterations
    verbose::Bool   # whether to show procedural information
    tol::Float64    # threshold for convergence
    lambda::Float64 # regularization parameter default 0

    function HierarchicalALS(;maxiter::Integer=100,
                             verbose::Bool=false,
                             tol::Real=1.0e-6,
                             lambda::Real=0.0)
        
        # check some constraints on the parameters.
        maxiter >= 1 || error("maxiter must be at least 1.")
        tol > 0 || error("tolerance must be positive.")
        # TODO: are there any constraints on lambda?
        #lambda >= || error("lambda must be nonnegative.")
    end
end

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


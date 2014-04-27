# Hierarchical Alternating Least Squares Algorithm for NMF.
# this file should be pulled into package JuliaStats/NMF.jl
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
export HierarchicalALS, solve!
immutable HierarchicalALS
    maxiter::Integer    # maximum number of iterations.
    tol::Float64        # threshold for convergence.
    lambda::Float64     # regularization parameter default 0.
    verbose::Bool       # whether to show procedural information.

    function HierarchicalALS(;maxiter::Integer=100,
                             tol::Real=1.0e-6,
                             lambda::Real=0.0,
                             verbose::Bool=false)
        
        # check some constraints on the parameters.
        maxiter >= 1 || error("maxiter must be at least 1.")
        tol > 0 || error("tolerance must be positive.")
        # TODO: are there any constraints on lambda?
        #lambda >= || error("lambda must be nonnegative.")
        new(maxiter, float(tol), float(lambda), verbose)
    end
end

function update_w!(W, AHt, HHt, k, epsilon)
    #to avoid divide by zero error.
    HHtDiag = diag(HHt)
    if !all(HHtDiag .> 0)
        throw(ArgumentError("HH' will prompt a division-by-zero"))
    end
    for x=1:k
        Wx = W[:,x] + (AHt[:,x] - W*HHt[:,x])/HHtDiag[x]
        Wx[Wx.<epsilon] = epsilon
        W[:,x] = Wx
    end
end

function update_h!(H, WtA, WtW, k, epsilon, lambda)
    #to avoid divide by zero error.
    WtWDiag=diag(WtW)
    if !all(WtWDiag .> 0)
        throw(ArgumentError("W'W will prompt a division-by-zero"))
    end
    for x=1:k
        #H[x,:] = H[x,:] + ((WtA[x,:] - WtW[x,:]*H) / WtWDiag[x])
        #H[x,:] = H[x,:] .- (lambda / WtWDiag[x])
        Hx = H[x,:] + ((WtA[x,:] - WtW[x,:]*H) / WtWDiag[x])
        Hx = Hx .- (lambda / WtWDiag[x])
        Hx[Hx .< epsilon] = epsilon
        H[x,:] = Hx
    end
    #H[H.<epsilon] = epsilon
end

function hals(A, W, H, k, tolerance, maxiter, lambda, verbose)
    # hals stands for Hierarchical Alternating Least Squares
    # solves A=WH
    # A = mxn matrix
    # W = mxk
    # H = kxn
    # k = rank of approximation.
    # implementation of algorithm2 from
    # http://www.bsp.brain.riken.jp/publications/2009/Cichocki-Phan-IEICE_col.pdf
    # All entries of A must be nonnegative.
    # A must not contain any rows with no zero elements.
    # We leave the user to decide how to ensure this.
    if !(all(sum(A, 1) .> 0) && all(sum(A, 2) .> 0))
        throw(ArgumentError("There is a zero row or column the algorithm will not converge."))
    end
    epsilon = eps(Float64)
    #prevError = vecnorm(A - W*H)
    oldW = copy(W)
    oldH = copy(H)
    prevError = vecnorm(W-oldW) + vecnorm(H-oldH)
    currError = prevError + 1
    currentIteration = 1
    errChange = zeros(1, maxiter)
    resids    = zeros(1, maxiter)
    converged = false
    while (!converged && currentIteration < maxiter)
        #update W
        AHt = A * H'
        HHt = H * H'
        copy!(oldW, W)
        update_w!(W, AHt, HHt, k, epsilon)

        #update H
        WtA = W' * A
        WtW = W' * W
        copy!(oldH, H)
        update_h!(H, WtA, WtW, k, epsilon, lambda)

        #check convergence
        if (currentIteration > 1)
            prevError = currError
        end
        errChange[currentIteration] = prevError
        #currError = vecnorm(A - W*H)
        currError = vecnorm(W-oldW) + vecnorm(H-oldH)
        if verbose
            resid = residual(A, W, H)
            resids[currentIteration] = resid
            if currentIteration > 1
                errChange[currentIteration] = resid - resids[currentIteration-1]
            else
                errChange[currentIteration] = resid
                @printf "iter\tchg W + H\tresid\tchgresid\n"
            end
            @printf("%d\t%f\t%f\t%f\n",
                    currentIteration,
                    currError,
                    resid,
                    errChange[currentIteration]
                    )
        end
        converged = abs(currError - prevError) < tolerance
        currentIteration = currentIteration + 1
    end
    result = NMF.Result(W, H, currentIteration, converged, currError)
    return result, resids, errChange
end

function solve!(alg::HierarchicalALS, A, W, H)
    m,n,k = NMF.nmf_checksize(A,W,H)
    result = hals(A, W, H, k, alg.tol, alg.maxiter, alg.lambda, alg.verbose)
    return result[1]
end

function randinit(m::Integer, n::Integer, k::Integer
                    ;normalize::Bool=false, zeroh::Bool=false)
   #X is m by n and we want a rank k factorization.
   W = rand(m, k)
   if normalize
       NMF.normalize1_cols!(W)
   end
   H = zeroh ? zeros(k, n) : rand(k, n)
   return (W, H)::(Matrix{Float64}, Matrix{Float64})
end

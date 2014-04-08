# Nonnegativity Constrained Least Squares with Multiple Righthand Sides 
#      using Block Principal Pivoting method
#
# This software solves the following problem: given A and B, find X such that
#              minimize || AX-B ||_F^2 where X>=0 elementwise.
#
# Reference:
#      Jingu Kim and Haesun Park, Toward Faster Nonnegative Matrix Factorization:
#           A New Algorithm and Comparisons,
#      In Proceedings of the 2008 Eighth IEEE International Conference on Data Mining 
#          (ICDM'08), 353-362, 2008
#
# Originally Written by Jingu Kim (jingu@cc.gatech.edu)
# Copyright 2008-2009 by Jingu Kim and Haesun Park, 
#                        School of Computational Science and Engineering,
#                        Georgia Institute of Technology
#
# Ported to Julialang by James Fairbanks (james.fairbanks@gatech.edu)
#             School of Computational Science and Engineering,
#             Georgia Institute of Technology
# Copyright 2014 by James Fairbanks, 
#                        School of Computational Science and Engineering,
#                        Georgia Institute of Technology
# Port has not been optimized for performance difference between Matlab and Julia.
#
# Last updated 2014-04-08 10:04
# Check updated code at http://www.github.com/jpfairbanks/Matfacgrf/
# Please send bug reports, comments, or questions to James Fairbanks.
# This code comes with no guarantee or warranty of any kind. 
# Note that this algorithm assumes that the
#      input matrix A has full column rank.
#
# Last modified 2014-04-08 10:00
#
# <Inputs>
#        A : input matrix (m x n) (by default), or A'*A (n x n) if isInputProd==1
#        B : input matrix (m x k) (by default), or A'*B (n x k) if isInputProd==1
#        isInputProd : (optional, default:0) if turned on, use (A'*A,A'*B) as input instead of (A,B)
#        init : (optional) initial value for X
# <Outputs>
#        X : the solution (n x k)
#        Y : A'*A*X - A'*B where X is the solution (n x k)
#        iter : number of iterations
#        success : 1 for success, 0 for failure.
#                  Failure could only happen on a numericall very ill-conditioned problem.

#require("solveNormalEqComb.jl")
export nnlsm_blockpivot
function nnlsm_blockpivot( A::Matrix, B::Matrix, isInputProd::Bool=false, init=false)
    if isInputProd
        AtA = A
        AtB = B
    else
        AtA = A'*A
        AtB = A'*B
    end
    
    n, k = size(AtB)
    MAX_ITER = n*5
    # set initial feasible solution
    X = zeros(n,k)
    if !init
        Y = - AtB
        PassiveSet = falses(n,k)
        iter = 0
    else
        PassiveSet = (init > 0)
        X, iter  = solveNormalEqComb(AtA,AtB,PassiveSet)
        Y = AtA * X - AtB
    end
    # parameters
    pbar = 3.0
    P = zeros(1,k)
    P[:] = pbar
    Ninf = zeros(1,k)
    Ninf[:] = n+1
    iter = 0

    NonOptSet = (Y .< 0.0) & ~PassiveSet
    InfeaSet  = (X .< 0.0) &  PassiveSet
    NotGood = sum(NonOptSet,1)+sum(InfeaSet,1) #relied on matlab summing over dim 1 by default
    NotOptCols = NotGood .> 0.0
    bigIter = 0
    success=1
    while(~isempty(find(NotOptCols)))
        bigIter = bigIter+1
        if ((MAX_ITER >0) && (bigIter > MAX_ITER))   # set max_iter for ill-conditioned (numerically unstable) case
            success = 0
            break
        end
        Cols1 = NotOptCols & (NotGood .< Ninf)
        Cols2 = NotOptCols & (NotGood .>= Ninf) & (P .>= 1)
        Cols3Ix = NotOptCols & ~Cols1 & ~Cols2
        Cols1Ints = find(Cols1)
        if ~isempty(Cols1Ints)
            P[Cols1Ints] = pbar
            Ninf[Cols1Ints] = NotGood[Cols1Ints]
            PassiveSet[find(NonOptSet & repmat(Cols1,n,1))] = true
            PassiveSet[find(InfeaSet & repmat(Cols1,n,1))] = false
        end
        Cols2Ints = find(Cols2)
        if ~isempty(Cols2Ints)
            P[Cols2Ints] = P[Cols2Ints]-1
            PassiveSet[NonOptSet & repmat(Cols2,n,1)] = true
            PassiveSet[InfeaSet & repmat(Cols2,n,1)] = false
        end
        Cols3IxInts = find(Cols3Ix)
        if ~isempty(Cols3IxInts)
            for i=1:length(Cols3IxInts)
                Ix = Cols3Ix[i]
                toChange = max(find( NonOptSet[:,Ix] | InfeaSet[:,Ix] )) #is | an or operator?
                if PassiveSet[toChange,Ix]
                    PassiveSet[toChange,Ix]=false
                else
                    PassiveSet[toChange,Ix]=true
                end
            end
        end
        NotOptMask = repmat(NotOptCols,n,1)
        NotOptColsInts = find(NotOptCols)
        X[:,NotOptCols], subiter = solveNormalEqComb(AtA,AtB[:,NotOptColsInts],PassiveSet[:,NotOptColsInts])
        iter = iter + subiter
        X[abs(X).<1e-12] = 0            # for numerical stability
        Y[:,NotOptCols] = AtA * X[:,NotOptColsInts] - AtB[:,NotOptColsInts]
        Y[abs(Y).<1e-12] = 0            # for numerical stability
        
        # check optimality
        NonOptSet = NotOptMask & (Y .< 0) & ~PassiveSet
        InfeaSet = NotOptMask & (X .< 0) & PassiveSet
        NotGood = sum(NonOptSet,1)+sum(InfeaSet,1)
        NotOptCols = NotGood .> 0
    end
    return X, Y, iter, success 
end

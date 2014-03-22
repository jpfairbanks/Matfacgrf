# utils.jl Utilities for graph analysis with NMF
# Author: James Fairbanks
# Date: 2014-03-17

#return the symmetric part of the graph/matrix
function symmetrize(A)
    return (A+A')/2
end

#normalize: a matrix along a dimension 1 means column stochastic.
function normalize(A, dim::Int)
    return A./sum(A,dim)
end

function residual(A,W,H)
    vecnorm(A - W*H)
end

# verify and return the sizes of the matrices
# A is m by n
# W is m by k
# H is k by n
# return m, k, n
function getsizes(A,W,H)
    m, n = size(A)
    k = size(W)[2]
    if k != size(H)[1] || n != size(H)[2]
         throw(DimensionMismatch("W, H are not compatible."))
    end
    return m, k, n
end

function rowresiduals(A,W,H,m,k,n)
    rowresids = zeros(n)
    for j = 1:n
        for i = 1:m
            x = 0.0
            #compute the predicted value
            for l = 1:k
                x += W[i,l] * H[l,j]
            end
            difference = (A[i,j] - x)
            rowresids[i] += difference^2
        end
    end
    return rowresids
end

function columnresiduals(A,W,H,m,k,n)
    colresids = zeros(n)
    for j = 1:n
        for i = 1:m
            x = 0.0
            #compute the predicted value
            for l = 1:k
                x += W[i,l] * H[l,j]
            end
            difference = (A[i,j] - x)
            colresids[j] += difference^2
        end
    end
    return colresids
end

function residual(A::SparseMatrixCSC, W, H, dimension::Integer)
    m, k, n = getsizes(A,W,H)
    if dimension == 2
        resids = columnresiduals(A, W, H, m, k, n)
    elseif dimension == 1
        resids = rowresiduals(A, W,H, m, k, n)
    else
        throw(ArgumentError("dimension must be 1 or 2."))
    end
end

function residual(A::SparseMatrixCSC, W, H)
    colresiduals = residual(A, W, H, 2)
    return sum(colresiduals)
end


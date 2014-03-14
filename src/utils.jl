
#return the symmetric part of the graph/matrix
function symmetrize(A)
    return (A+A')/2
end

#normalize: a matrix along a dimension 1 means column stochastic.
function normalize(A, dim::Int)
    return A./sum(A,dim)
end


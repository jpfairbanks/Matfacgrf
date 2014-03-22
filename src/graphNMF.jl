using MLBase
using NMF

# statc matrix factorization on entire graph with rank k
# Returns the matrix representing the graph and the NMF.Result.
function graphNMF(alg::HierarchicalALS, AdjMat, k::Int)
    AdjMat += speye(size(AdjMat)[1])
    S = symmetrize(AdjMat)
    X = S
    #X = normalize(S, 1)
    W, H = randinit(size(X)[1], size(X)[2], k)
    #nmfresult = hals(AdjMat, W, H, k, tolerance, 50,0)
    nmfresult = solve!(alg, X, W, H)
    return X, nmfresult
end

# make a closure that will accept updates to 
# a sparse matrix and return the updated factors
# as an NMF.Result
function NMFClosure(alg::HierarchicalALS, maxVertices::Int, rank::Int)
    # initialize
    AdjMat = speye(maxVertices,maxVertices)
    #do random initialization once to improve continuity
    W, H = randinit(maxVertices, maxVertices, rank)
    function batchHandle(batch)
        AdjMat += batch
        #result = hals(AdjMat, W, H, rank, tolerance, 50,0)
        result = solve!(alg, AdjMat, W, H)
        return result
    end
    return batchHandle
end

# Run nmf of a graph and classify the vertices into their clusters.
# Returns the labels and the number of vertices in each class.
function nmfclassify(alg::HierarchicalALS, AdjMat, k)
    A, result = graphNMF(alg, AdjMat, k)
    H = result.H
    labels = MLBase.classify(H)
    counts = hist(labels)
    return labels, counts
end


#compute the residuals from a rank k approximation.
function nmfresiduals(alg::HierarchicalALS, AdjMat, k::Integer)
    X, result = graphNMF(alg, AdjMat, k)
    W = result.W
    H = result.H
    return residual(X, W, H, 2)
end

# TODO: Implement a function that compares the effects
# of normalization on the outputs of the NMF on graphs.
# Does normalizing by degree help or hurt us for resolving
# clusters or patterns?
function compareNormalizations()
    error("Not implemented yet")
end

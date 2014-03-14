# hospital.jl : Loads the hospital contact graph and prints the adjacency matrix
# Author: James Fairbanks
# Date:   2014-03-12

include("../src/utils.jl")
include("../src/stream_load.jl")
include("../src/hals.jl")
using MLBase
using Gadfly

const tolerance = 0.00001
dataset = FileParams(
    "data/hospital_edges.csv",
    100,
    75)

#statc matrix factorization on entire graph with rank k
function graphNMF(k::Int)
    AdjMat = readgraph(dataset)
    AdjMat += speye(size(AdjMat)[1])
    S = symmetrize(AdjMat)
    X = normalize(S, 1)
    W, H = randinit(size(X)[1], size(X)[2], k)
    W, H = hals(AdjMat,W,H,k, tolerance, 50,0)
    return X, W, H
end

#histogram the residuals from a rank k approximation.
function histresiduals(k::Int)
    A, W, H = graphNMF(k)
    plt = plot(x=residuals(A,W,H), Geom.histogram)
    draw(SVG("hist.svg", 12cm,12cm), plt)
end


function NMFClosure(maxVertices::Int, rank::Int)
    # initialize
    Adjmat = speye(maxVertices,maxVertices)
    X = full(Adjmat)
    #do random initialization once to improve continuity
    W, H = NMF.randinit(maxVertices, maxVertices, rank)
    function batchHandle(batch)
        Adjmat += batch
        W, H = hals(Adjmat,W,H,rank, tolerance, 50,0)
        return H
    end
    return batchHandle
end


function hospital_demo()
    @time begin
        handler = showClosure(dataset.maxVertices)
        processBatches(dataset, handler)
    end
end

function plotverts(H)
    p = plot(x=H[1,:], y=H[2,:])
    print("Writing to file hospital.svg\n")
    draw(SVG("hospital.svg", 6inch,3inch), p)
end

#dynamic factorization operating on accumulating graph.
# Returns an H for each timestep
function dynamic_graphNMF()
    k = 4
    #H = zeros(k,dataset.maxVertices)
    i = 0
    handler = NMFClosure(dataset.maxVertices, k)
    df = readtable(dataset.file)
    time = cell(iceil(size(df)[1]/dataset.batchsize))
    batches = @task yieldBatchMats(df, dataset.batchsize, dataset.maxVertices)
    for M in batches
        i += 1
        H = handler(M)
        time[i] = H./sum(H,1)
    end
    return time
end

function hospital_classify()
    A, W, H = graphNMF(4)
    labels = MLBase.classify(H)
    counts = hist(labels)
end


#static version operating on the full graph. Makes 2D plot.
function hospital_plot_vertices()
    A, W, H = graphNMF(2)
    plotverts(H)
end

# hospital.jl : Loads the hospital contact graph and prints the adjacency matrix
# Author: James Fairbanks
# Date:   2014-03-12

using Matfacgrf
using MLBase
using Winston
using DataFrames

import Matfacgrf.FileParams
import Matfacgrf.HierarchicalALS
import Matfacgrf.readgraph
import Matfacgrf.graphNMF
import Matfacgrf.nmfresiduals
import Matfacgrf.nmfclassify
import Matfacgrf.NMFClosure
import Matfacgrf.yieldBatchMats

const tolerance = 0.00001

# static function to assign a class label to each vertex
function classifyvertices(alg, AdjMat, k::Integer)
    labels, counts = nmfclassify(alg, AdjMat, k)
end

##############################
# Plotting Code              #
##############################

#histogram the residuals from a rank k approximation.
function histresiduals(filename, k::Int)
    resids = nmfresiduals(alg, AdjMat, k)
    plt = plothist(resids)
    info("Drawing distribution of residuals to file: $filename\n")
    file(filename)
end

function plotcolumns(filename, H,)
    p = scatter(H[1,:], H[2,:])
    info("Drawing vertex embedding to file: $filename\n")
    file(filename)
end

#static version operating on the full graph. Makes 2D plot.
function plotvertices(filename, alg, AdjMat)
    X, result = graphNMF(alg, AdjMat, 2)
    H = result.H
    plotcolumns(filename, H)
end


################################
#    Streaming Code            #
################################
function batchcat()
    @time begin
        handler = showClosure(dataset.maxVertices)
        processBatches(dataset, handler)
    end
end
#dynamic factorization operating on accumulating graph.
# Returns an H for each timestep
function dynamic_graphNMF(dataset::FileParams, k::Integer)
    #H = zeros(k,dataset.maxVertices)
    i = 0
    handler = NMFClosure(alg, dataset.maxVertices, k)
    df = DataFrames.readtable(dataset.file)
    time = cell(iceil(size(df)[1]/dataset.batchsize))
    batches = @task yieldBatchMats(df, dataset.batchsize, dataset.maxVertices)
    for M in batches
        i += 1
        H = handler(M).H
        time[i] = H./sum(H,1)
    end
    return time
end



dataset = FileParams(
    "data/hospital_edges.csv",
    3000,
    75,
    2,
    3)

info("Reading graph from $dataset.")
AdjMat = readgraph(dataset)
alg = HierarchicalALS()

function testHospital(k::Integer)
    info("A static embedding of the vertices based on NMF.")
    vertexplot("hospscatter.svg", alg, AdjMat)
    info("You can find outliers based on the residuals")
    histresiduals("histogram.svg", k)
    info("Classifying vertices into $k groups.")
    @show vertexclassify(alg, AdjMat, k)
    info("We can update the embedding at each batch.")
    locations = dynamic_graphNMF(dataset, k)
    info("Test finished")
end

testHospital(6)

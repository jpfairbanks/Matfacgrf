# visualize.jl : Loads the graph and make plots.
# Author: James Fairbanks
# Date:   2014-03-12

using Matfacgrf
using MLBase
using Winston
using DataFrames
using ArgParse

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
function histresiduals(filename, alg, AdjMat, k::Int)
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
function dynamic_graphNMF(dataset::FileParams, alg, k::Integer)
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

function getsettings()
    args = ArgParseSettings()
    @add_arg_table args begin
        "--file", "-f"
            help = "file containing edges"
            arg_type = String
        "--batchsize"
            help = "number of edges per batch"
            arg_type = Int
        "--max_vertices"
            help = "largest vertex ID that will be encountered"
            arg_type = Int
        "--rank", "-k"
            help = "rank of the factorization"
            arg_type = Int
    end
    settings = parse_args(args)
    return settings
end

function main(args)
    dataset = FileParams(
        args["file"],
        args["batchsize"],
        args["max_vertices"],
        2,
        3)
    k = args["rank"]
    @show dataset
    info("Reading graph from $dataset.")
    AdjMat = readgraph(dataset)
    alg = HierarchicalALS()
    info("A static embedding of the vertices based on NMF.")
    plotvertices("vertexscatter.svg", alg, AdjMat)
    info("You can find outliers based on the residuals")
    histresiduals("residualhistogram.svg", alg, AdjMat, k)
    info("Classifying vertices into $k groups.")
    @show classifyvertices(alg, AdjMat, k)
    info("We can update the embedding at each batch.")
    locations = dynamic_graphNMF(dataset, alg, k)
end

settings = getsettings()
main(settings)

# stream_load.jl : Loads a graph in batches from a file using DataStreams
# Author: James Fairbanks
# Date:   2014-03-12
using DataFrames

function showClosure(maxVertices)
    Adjmat = spzeros(maxVertices,maxVertices)
    function batchHandle(batch)
        @show batch
        Adjmat += batch
        @show Adjmat
    end
    return batchHandle
end

function yieldBatchMats(datastream, maxVertices)
    for df in datastream
        #since the first column is the row number
        #convert from python's 0 indexing to 1 indexing
        batchsrc = convert(Array, df[2] + 1)
        batchdest = convert(Array, df[3] + 1)
        batchMat = sparse(batchsrc, batchdest, 1, maxVertices, maxVertices)
        produce(batchMat)
    end
end
function processBatches(filename, batchsize, maxVertices, handle)
    datastream = readstream(edgefile, nrows=batchsize)
    batchMats = @task yieldBatchMats(datastream, maxVertices)
    for M in batchMats
        handle(M)
    end
end

function hospital_main()
    edgefile = "data/hospital_edges.csv"
    batchsize = 100
    maxV = 75
    handler = showClosure(maxV)
    processBatches(edgefile, batchsize,maxV, handler)
end

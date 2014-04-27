# stream_load.jl : Loads a graph in batches from a file using DataStreams
# Author: James Fairbanks
# Date:   2014-03-12
using DataFrames

immutable FileParams
    file
    batchsize::Integer
    maxVertices::Integer
    srccol::Integer
    dstcol::Integer
end

function readgraph(fp::FileParams)
    df = readtable(fp.file)
    #since the first column is the row number
    #convert from python's 0 indexing to 1 indexing
    src  = df.columns[fp.srccol] .+ 1
    dest = df.columns[fp.dstcol] .+ 1
    m = maximum([maximum(src), maximum(dest)])
    AdjMat = sparse(src, dest, [1.0 for i in src], m,m)
    return AdjMat
end

function yieldBatchMats(dataframe, batchsize, maxVertices)
    #since the first column is the row number
    #convert from python's 0 indexing to 1 indexing
    src = convert(Array, dataframe[1] .+ 1)
    dst = convert(Array, dataframe[2] .+ 1)
    for i = 1:batchsize:size(dataframe)[1]
        starti = i
        endi   = min(i+batchsize, length(src))
        batchsrc = src[starti:endi]
        batchdst = dst[starti:endi]
        batchMat = sparse(batchsrc, batchdst, 1, maxVertices, maxVertices)
        produce(batchMat)
    end
end

function yieldBatchMats(fp::FileParams)
    df = readtable(fp.file)
    stream = @task yieldBatchMats(df, fp.batchsize, fp.maxVertices)
    return stream
end

function showClosure(maxVertices)
    Adjmat = spzeros(maxVertices,maxVertices)
    function batchHandle(batch)
        @show batch
        Adjmat += batch
        #@show Adjmat
    end
    return batchHandle
end

function processBatches(filep::FileParams, handle)
    datastream = DataFrames.readtable(filep.file)
    batchMats = @task yieldBatchMats(datastream, filep.batchsize, filep.maxVertices)
    for M in batchMats
        handle(M)
    end
end

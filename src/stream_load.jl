using DataFrames

function produceMonad(maxVertices)
    Adjmat = spzeros(maxVertices,maxVertices)
    function batchHandle(batch)
        @show batch
        Adjmat += batch
        @show Adjmat
    end
    return batchHandle
end

function processBatches(filename, batchsize, maxVertices, batchMonad)
    datastream = readstream(edgefile, nrows=batchsize)
    for df in datastream
        #since the first column is the row number
        #convert from python's 0 indexing to 1 indexing
        batchsrc = convert(Array, df[2] + 1)
        batchdest = convert(Array, df[3] + 1)
        batchMat = sparse(batchsrc, batchdest, 1, maxVertices, maxVertices)
        batchMonad(batchMat)
    end
end

edgefile = "data/hospital_edges.csv"
batchsize = 100
maxV = 75
handler = produceMonad(maxV)
processBatches(edgefile, batchsize,maxV, handler)

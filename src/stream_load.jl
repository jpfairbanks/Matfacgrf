using DataFrames

edgefile = "data/hospital_edges.csv"
batchsize = 100
maxV = 75
datastream = readstream(edgefile, nrows=batchsize)
Adjmat = spzeros(maxV,maxV)
for df in datastream
    #since the first column is the row number
    batchsrc = convert(Array, df[2] + 1)
    batchdest = convert(Array, df[3] + 1)
    batchMat = sparse(batchsrc, batchdest, 1, maxV, maxV)
    @show batchMat
    Adjmat += batchMat
    @show Adjmat
end

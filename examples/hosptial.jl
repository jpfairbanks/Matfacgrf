# hospital.jl : Loads the hospital contact graph and prints the adjacency matrix
# Author: James Fairbanks
# Date:   2014-03-12

include("../src/stream_load.jl")

function hospital_main()
    edgefile = "data/hospital_edges.csv"
    batchsize = 100
    maxV = 75
    handler = showClosure(maxV)
    processBatches(edgefile, batchsize,maxV, handler)
end

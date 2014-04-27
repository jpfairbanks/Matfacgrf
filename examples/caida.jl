using Matfacgrf
#using Winston
#using DataFrames
using ArgParse
using MAT

import Matfacgrf.FileParams
import Matfacgrf.HierarchicalALS
import Matfacgrf.readgraph
import Matfacgrf.graphNMF

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

function report_sparsity(W)
    @show size(W)[1]
    @printf "W.nfilled, W.nnz, density : %d, %d, %f\n" nfilled(W) countnz(W) countnz(W)/length(W)
    @show sum(W, 1)
    @show sum([W .!= 0], 1)
end

function loaddata(args)
    dataset = FileParams(
        args["file"],
        args["batchsize"],
        args["max_vertices"],
        2,
        3)
    AdjMat = readgraph(dataset, dataset.batchsize)
    k = args["rank"]
    println("loaded graph")
    return dataset, AdjMat, k
end

function main(args)
    # input time
    dataset, AdjMat, k = @time loaddata(args)
    alg = HierarchicalALS()
    
    #compile
    graphNMF(alg, speye(5,5), k)

    # do the work
    X, result = @time graphNMF(alg, AdjMat, k)
    #writedlm("H.mat", result.H, " ")
    W = result.W

    # normalize columns of W
    W = W ./ sum(W,1)

    # just drop the really small entries in W
    W[W.<1e-10] = 0
    report_sparsity(W)

    #output
    writedlm("W.mat", W, " ")
    writedlm("A.mat", readgraph(dataset), " ")
end


# matlabIO: take and NMF result and dump it to a pair of mat files.
function matlabIO(result)
    file = matopen("W.mat", "w")
    write(file, "W", scaledW)
    close(file)
    file = matopen("H.mat", "w")
    write(file, "H",result.H)
    close(file)
end

settings = getsettings()
print("got settings")
print(settings)
main(settings)

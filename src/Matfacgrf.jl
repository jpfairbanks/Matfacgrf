module Matfacgrf
    include("stream_load.jl")
    include("utils.jl")
    include("hals.jl")
    include("graphNMF.jl")
    include("solveNormalEqComb.jl")
    include("nnlsm_blockpivot.jl")
end

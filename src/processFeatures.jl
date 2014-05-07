using Autoreload
arequire("derivedFeatures.jl")

function main(args)
    extractors = [id! paddedDiffsqr! movAvg(0.5) knuthVar!]
    @show extractors
    for filename in args
        @time F = readdlm(filename)
        print("Read: $filename\n")
        M = derivedFeatures(F, extractors)
        @time writedlm("$filename.feat", M)
        #@show M
    end
end

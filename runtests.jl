using Matfacgrf

tests = ["utils",
         "hals",
         "recover",
         "nnls"]

println("Running tests:")

for t in tests
    tfile = joinpath("test", "$(t).jl")
    println(" * $(tfile) ...")
    include(tfile)
end

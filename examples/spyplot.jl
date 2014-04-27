using Gadfly

filename = "W.mat"
ext = ".pdf"
fp = open(filename, "r")
W = readdlm(fp, ' ', '\n')
#Winston.spy(W)
almostnonzero(x) = x > 1e-5
linindex = find(almostnonzero, W)
rowid = div(linindex, size(W)[1]) .+ 1
colid = rem(linindex, size(W)[1]) .+ 1
val = W[linindex]
draw(SVG("spy.svg", 6inch, 3inch), Gadfly.spy(sparse(rowid, colid, val)))
#Winston.file(joinpath(filename, ext))

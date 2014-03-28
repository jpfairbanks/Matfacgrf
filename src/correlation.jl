require("src/stream_load.jl")
require("src/utils.jl")
require("src/hals.jl")

function pagerank(G::Array{Float64,2}, tele::Float64, initial::Array{Float64,1})
	# computes the pagerank vector of a dense Graph G with teleportation constant tele
	# and an initial guess of initial. If initial does not have any projection onto 
	# the pagerank vector it will not be correct.
	A = G + (tele*ones(size(G)))
	A = normalize(A, 1)
	p = initial ./sum(initial)
	for i in 1:100
		p = A*p
	end
	
	evd = eig(A)
	pageranks = evd[2][:,1]
	#@show pageranks
	@assert all(sign(p).==1) || all(sign(p).==-1) "pageranks are not all the same sign!"
	@assert all(sign(pageranks).==1) || all(sign(pageranks).==-1) "pageranks are not all the same sign!"
	#x = pageranks ./ sum(pageranks)
	return p
end

function maindegree(fp::FileParams)
	numbatches = 11
	G = eye(fp.maxVertices, fp.maxVertices)
	#pr = rand(fp.maxVertices)
	stream = yieldBatchMats(fp)
	X = zeros(fp.maxVertices, numbatches)
	deg = sum(G, 2)
	for i in 1:numbatches
		batch = consume(stream)
		G += (batch+batch')/2
		sum!(deg, G)
		X[:, i] = deg
	end
	return X
end

function mainpagerank(fp::FileParams)
	tele = 0.05 #pagerank teleportation constant
	G = eye(fp.maxVertices, fp.maxVertices)
	#pr = rand(fp.maxVertices)
	stream = yieldBatchMats(fp)
	numbatches = 11
	X = zeros(fp.maxVertices, numbatches)
	for i = 1:numbatches
		batch = consume(stream)
		G += (batch+batch')/2
		X[:, i] = pagerank(G, tele, rand(fp.maxVertices))
	end
	return X
end

function mainhals(fp::FileParams)
	G = eye(fp.maxVertices, fp.maxVertices)
	stream = yieldBatchMats(fp)
	numbatches = 11
	X = zeros(fp.maxVertices, numbatches)
	alg = HierarchicalALS()
	for i = 1:numbatches
		batch = consume(stream)
		G += (batch+batch')/2
		W, H = randinit(fp.maxVertices, fp.maxVertices, 4)
		res = solve!(alg, G, W, H,)
		H = res.H
		nH = H./sum(H,2)
		X[:, i] = [vecnorm(nH[:,i]) for i in 1:fp.maxVertices]
	end
	return X
end

function test()
	hospital = FileParams(
		"./data/hospital_edges.csv",
		3000,
		75,
		1,2
	)
	#println(string(mainpagerank(hospital)))
	#vtxfunc(h) = vecnorm2(h)
	degrees_over_time= maindegree(hospital)
	@show cor(degrees_over_time)
	pageranks_over_time = mainpagerank(hospital)
	@show cor(pageranks_over_time)
	hals_uncertainty_over_time = mainhals(hospital)
	@show cor(hals_uncertainty_over_time)
end

using Base.Test 
using Matfacgrf

function testnnls(m::Integer,n::Integer,k::Integer, tol::Float64)
	W = rand(m,k)
	H = rand(k,n)
	A = W*H
	Result = Matfacgrf.nnlsm_blockpivot(W,A)
	Hhat = Result[1]
	@test norm(H-Hhat) <= tol
end

testnnls(8,10,4, 1e-10)
testnnls(16,10,4, 1e-10)
testnnls(8,10,8, 1e-10)
@test_throws testnnls(4,4,5,1e-10)


using Matfacgrf
#import Matfacgrf.randinit
#import Matfacgrf.HierarchicalALS
#import Matfacgrf.solve!

using Base.Test
function testHALStype()
    @test_throws HierarchicalALS(maxiter = 0)
    @test_throws HierarchicalALS(maxiter= -1)
    @test_throws HierarchicalALS(tol= -10e-6)
    @test_throws HierarchicalALS(verbose=1)
    @test isa(HierarchicalALS(lambda= -1), HierarchicalALS)
    @test isa(HierarchicalALS(verbose=true), HierarchicalALS)
end

function testHALS()
    m, n, k = 10, 15, 5
    W, H = Matfacgrf.randinit(m, n, k)
    A = W * H + rand(m,n)/5
    alg = HierarchicalALS(maxiter=5,
                           tol=10e-6,
                           lambda=0.0,
                           verbose=true)
    @test isa(solve!(alg, A, W, H), NMF.Result)
    alg = HierarchicalALS(maxiter=500,
                           tol=10e-6,
                           lambda=0.0,
                           verbose=false)
    res = solve!(alg, A, W, H)
    @test isa(res, NMF.Result)
    @show res.niters, res.converged, res.objvalue
end

testHALStype()
testHALS()

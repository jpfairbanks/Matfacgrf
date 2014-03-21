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
    W, H = randinit(m, n, k)
    A = W * H + rand(m,n)/5
    alg = HierarchicalALS(maxiter=100,
                           tol=10e-6,
                           lambda=0.0,
                           verbose=false)
    @test isa(solve!(alg, A, W, H), NMF.Result)
end

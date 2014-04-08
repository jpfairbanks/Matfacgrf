#adaptation of Base.sortrows builtin to return the permutation like matlab does
import Base.Order.Lexicographic
function sortrowsperm(A::AbstractMatrix; kws...)
    c = 1:size(A,2)
    rows = [ sub(A,i,c) for i=1:size(A,1) ]
    p = sortperm(rows; kws..., order=Lexicographic)
    return A[p,:], p
end

function solveNormalEqComb( AtA, AtB, PassSet )
# Solve normal equations using combinatorial grouping.
# Although this function was originally adopted from the code of
# "M. H. Van Benthem and M. R. Keenan, J. Chemometrics 2004; 18: 441-450",
# important modifications were made to fix bugs.
#
# Modified by Jingu Kim (jingu@cc.gatech.edu)
#             School of Computational Science and Engineering,
#             Georgia Institute of Technology
#
# Ported to Julialang by James Fairbanks (james.fairbanks@gatech.edu)
#             School of Computational Science and Engineering,
#             Georgia Institute of Technology
# Last updated 2014-04-08 10:04
    iter = 0;
    if isempty(PassSet) || all(PassSet)
        Z = AtA\AtB
        iter = iter + 1
    else
        Z = zeros(size(AtB))
        n, k1 = size(PassSet)
        PassSetInt = find(PassSet)
        ## Fixed on Aug-12-2009
        if k1==1
            Z[PassSetInt]=AtA[PassSetInt,PassSetInt]\AtB[PassSetInt]
        else
            ## Fixed on Aug-12-2009
            # The following bug was identified by investigating a bug report by Hanseung Lee.
            sortedPassSet, sortIx = sortrowsperm(PassSet')
            breaks = any(bool(diff(sortedPassSet,1)'))
            breakIx = [0 find(breaks) k1]
            # codedPassSet = 2.^(n-1:-1:0)*PassSet
            # [sortedPassSet,sortIx] = sort(codedPassSet)
            # breaks = diff(sortedPassSet)
            # breakIx = [0 find(breaks) k1]

            for k=1:length(breakIx)-1
                cols = sortIx[breakIx[k]+1:breakIx[k+1]]
                vars = PassSet[:,sortIx[breakIx[k]+1]]
                Z[vars,cols] = AtA[vars,vars]\AtB[vars,cols]
                iter = iter + 1
            end
        end
    end
    return Z, iter
end

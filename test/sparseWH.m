function [W,H] = sparseWH(m,k,n,densityW, densityH)
    %make a random sparse W,H pair for testing.
    W = sprand(m, k, densityW);
    H = sprand(k, n, densityH);
end

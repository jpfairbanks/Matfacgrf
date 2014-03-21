function [W,H] = denseWH(m,k,n)
    %make a W,H pair that is random for testing purposes.
    W = rand(m, k);
    H = rand(k, n);
end

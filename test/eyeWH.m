function [W, H] = eyeWH(m,k,n)
%make a W,H pair that is simple for testing purposes.
    W = eye(m,k);
    H = eye(k,n);
end

function [beta, x] = inviter(A, s, numiter)
% INVITER   Shifted inverse iteration for the closest eigenvalue.
% Input:
%   A         square matrix
%   s         value close to targeted eigenvalue (complex scalar)
%   numiter   number of iterations
% Output: 
%   beta      sequence of eigenvalue approximations (vector)
%   x         final eigenvector approximation

    n = length(A);
    x = randn(n, 1);
    x = x / norm(x, inf);
    B = A - s * eye(n);
    [L, U] = lu(B);
    beta = zeros(numiter, 1);
    for k = 1:numiter
        y = U \ (L \ x);
        beta(k) = (1 / dot(x, y)) + s;
        x = y / norm(y);
    end
end

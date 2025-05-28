function [beta, x] = poweriter(A, numiter)
% POWERITER   Power iteration for the dominant eigenvalue.
% Input:
%   A         square matrix
%   numiter   number of iterations
% Output: 
%   beta      sequence of eigenvalue approximations (vector)
%   x         final eigenvector approximation

    n = length(A);
    x = randn(n, 1);
    x = x / norm(x);
    beta = zeros(numiter, 1);
    for k = 1:numiter
        y = A * x;
        beta(k) = x' * y;
        x = y / norm(y);
    end
end
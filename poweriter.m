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
x = x / norm(x, inf);
for k = 1:numiter
    y = A*x;
    [normy, m] = max(abs(y));
    beta(k) = y(m) / x(m);
    x = y / y(m);
end 
function x = newtonsys(f,x1)
% NEWTONSYS   Newton's method for a system of equations.
% Input:
%   f        function that computes residual and Jacobian matrix
%   x1       initial root approximation (n-vector)
% Output       
%   x        array of approximations (one per column, last is best)

    % Stopping parameters.
    funtol = 1000 * eps;    % stop for small || f(x) ||
    xtol = 1000 * eps;      % stop for small || x change ||
    maxiter = 40;           % stop after this many iterations

    x = x1(:);  
    [y, J] = f(x1);
    dx = Inf;
    k = 1;

    while (norm(dx) > xtol) && (norm(y) > funtol) && (k < maxiter)
        dx = -(J \ y);    % Newton step
        x(:, k+1) = x(:, k) + dx;

        k = k+1;
        [y, J] = f(x(:, k));
    end

    if k==maxiter
        warning('Maximum number of iterations reached.')
    end
end
function x = newton(f, df_dx, x1)
% NEWTON   Newton's method for a scalar equation.
% Input:
%   f        objective function 
%   df_dx    derivative function
%   x1       initial root approximation
% Output       
%   x        vector of root approximations (last one is best)

    % Stopping parameters.
    funtol = 100 * eps;    % stop for small f(x)
    xtol = 100 * eps;      % stop for small x change
    maxiter = 40;          % stop after this many iterations

    x = x1;  
    y = f(x1);
    dx = Inf;   % for initial pass below
    k = 1;

    while (abs(dx) > xtol) && (abs(y) > funtol) && (k < maxiter)
        dydx = df_dx(x(k));
        dx = -y / dydx;    % Newton step
        x(k+1) = x(k) + dx;

        k = k+1;
        y = f(x(k));
    end

    if k==maxiter
      warning('Maximum number of iterations reached.')
    end
end
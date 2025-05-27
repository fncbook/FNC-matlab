function x = secant(f, x1, x2)
% SECANT   Secant method for a scalar equation.
% Input:
%   f        objective function 
%   x1, x2   initial root approximations
% Output       
%   x        vector of root approximations (last is best)

    % Stopping parameters.
    funtol = 100 * eps;    % stop for small f(x)
    xtol = 100 * eps;      % stop for small x change
    maxiter = 40;          % stop after this many iterations

    x = [x1 x2];
    y1 = f(x1);  y2 = f(x2);
    k = 2;       % iteration counter
    dx = Inf;    % for initial pass below

    while (abs(dx) > xtol) && (abs(y2) > funtol) && (k < maxiter)
        dx = -y2 * (x(k) - x(k-1)) / (y2 - y1);    % secant step
        x(k+1) = x(k) + dx;
        
        k = k+1;
        y1 = y2;    % current f-value becomes the old one next time
        y2 = f(x(k));
    end

    if k==maxiter
        warning('Maximum number of iterations reached.')
    end
end

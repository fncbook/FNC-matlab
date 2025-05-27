function [x, u, du_dx] = shoot(phi, a, b, ga, gb, init, tol)
%SHOOT   Shooting method for a two-point boundary-value problem.
% Input:
%   phi      defines u'' = phi(x, u, u') (function)
%   a, b     endpoints of the domain (scalars)
%   ga       residual boundary function of u(a), u'(a) 
%   gb       residual boundary function of u(b), u'(b) 
%   init     initial guess for u(a) and u'(a) (column vector)
%   tol      error tolerance (scalar)
% Output:
%   x        nodes in x (length n+1)
%   u        values of u(x)  (length n+1)
%   du_dx    values of u'(x) (length n+1)

    % To be solved by the IVP solver
    function f = shootivp(x, y)
      f = [ y(2); phi(x, y(1), y(2)) ];
    end

    % To be solved by levenberg
    function residual = objective(s)
      [x, y] = rk23(@shootivp, [a, b], s, tol / 10);
      ya = y(1, :);
      yb = y(end, :);
      residual = [ga(ya(1), ya(2)); gb(yb(1), yb(2))];
    end

    y = [];    % shared variable
    s = levenberg(@objective, init, tol);

    % Don't need to solve the IVP again. It was done within the
    % objective function already.
    u = y(:, 1);            % solution     
    du_dx = y(:, 2);        % derivative
end
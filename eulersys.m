function [t,u] = eulersys(du_dt, tspan, u0, n)
% EULERSYS   Euler's method for a first-order IVP system.
% Input:
%   du_dt   defines f in u'(t) = f(t, u) 
%   tspan   endpoints of time interval (2-vector)
%   u0      initial value (m-vector)
%   n       number of time steps (integer)
% Output:
%   t       selected nodes (vector, length n+1)
%   u       solution values (array, n+1 by m)

    % Define the time discretization.
    a = tspan(1);  b = tspan(2);
    h = (b - a) / n;
    t = a + (0:n)' * h;

    % Initialize solution array.
    u = zeros(length(u0), n+1);
    u(:, 1) = u0(:);

    % Time stepping.
    for i = 1:n
      u(:, i+1) = u(:, i) + h * du_dt(t(i), u(:, i));
    end

    u = u.';    % conform to MATLAB output convention
end
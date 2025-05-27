function [t, u] = ab4(du_dt, tspan, u0, n)
%AB4     4th-order Adams-Bashforth formula for an IVP.
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

    % Constants in the AB4 method.
    k = 4;  
    sigma = [55; -59; 37; -9] / 24;  

    % Find starting values by RK4.
    u = zeros(length(u0), n+1);
    [ts, us] = rk4(du_dt,[a, a + (k-1)*h], u0, k-1);
    u(:, 1:k) = us(1:k, :).';

    % Compute history of u' values, from oldest to newest.
    f = zeros(length(u0), k);
    for i = 1:k-1
      f(:, k-i) = du_dt(t(i), u(:, i));
    end

    % Time stepping.
    for i = k:n
      f = [du_dt(t(i), u(:, i)), f(:, 1:k-1)];    % new value of du/dt
      u(:, i+1) = u(:, i) + h * (f * sigma);      % advance one step
    end

    u = u.';   % conform to MATLAB output convention
end
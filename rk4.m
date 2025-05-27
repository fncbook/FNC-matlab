function [t, u] = rk4(du_dt, tspan, u0, n)
% RK4    Fourth-order Runge-Kutta for an IVP.
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
      k1 = h * du_dt( t(i),       u(:, i)       );
      k2 = h * du_dt( t(i) + h/2, u(:, i) + k1/2);
      k3 = h * du_dt( t(i) + h/2, u(:, i) + k2/2);
      k4 = h * du_dt( t(i) + h,   u(:, i) + k3  );
      u(:, i+1) = u(:, i) + (k1 + 2*(k2 + k3) + k4) / 6;
    end

    u = u.';    % conform to MATLAB output convention
end
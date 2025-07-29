function [t, u] = am2(du_dt, tspan, u0, n)
    % AM2    2nd-order Adams-Moulton (trapezoid) formula for an IVP.
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

    u = zeros(length(u0), n+1);
    u(:, 1) = u0(:);

    % This function defines the rootfinding problem at each step.
    function F = trapzero(z)
        F = z - h/2 * du_dt(t(i+1), z) - known;
    end

    % Time stepping.
    for i = 1:n
        % Data that does not depend on the new value.
        known = u(:,i) + h/2 * du_dt(t(i), u(:, i));
        % Find a root for the new value. 
        unew = levenberg(@trapzero, known);
        u(:, i+1) = unew(:, end);
    end

    u = u.';   % conform to MATLAB output convention
end  % main function
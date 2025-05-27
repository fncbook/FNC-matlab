function p = plinterp(t, y)
% PLINTERP   Piecewise linear interpolation.
% Input:
%   t     interpolation nodes (vector, length n+1)
%   y     interpolation values (vector, length n+1)
% Output:
%   p     piecewise linear interpolant (function)

    n = length(t) - 1;
    H = {};
    for k = 0:n
        H{k+1} = hatfun(t, k);
    end
    p = @evaluate;

    % This function evaluates p when called.
    function f = evaluate(x)
        f = 0;
        for i = 0:n
            f = f + y(i+1) * H{i+1}(x);
        end
    end
end
function H = hatfun(t, k)
% HATFUN   Hat function/piecewise linear basis function.
% Input: 
%   t      interpolation nodes (vector, length n+1)
%   k      node index (integer, in 0,...,n)
% Output:
%   H      kth hat function (function)

    n = length(t) - 1;
    H = @evaluate;

    function y = evaluate(x)
        y = zeros(size(x));
        for i = 1:numel(x)
            if (k > 0) && (t(k) <= x(i)) && (x(i) <= t(k+1))
                y(i) = (x(i) - t(k)) / (t(k+1) - t(k));
            elseif (k < n) && (t(k+1) <= x(i)) && (x(i) <= t(k+2))
                y(i) = (t(k+2) - x(i)) / (t(k+2) - t(k+1));
            end
        end
    end
end
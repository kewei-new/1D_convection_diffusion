function orders = compute_orders_from_errors(errors)
%COMPUTE_ORDERS_FROM_ERRORS Compute convergence orders from error vector.
%
% Input:
%   errors : column vector of errors
%
% Output:
%   orders : column vector, first entry is NaN

    errors = errors(:);
    n = length(errors);

    orders = nan(n,1);

    for k = 2:n
        if errors(k) > 0 && errors(k-1) > 0
            orders(k) = log2(errors(k-1) / errors(k));
        else
            orders(k) = NaN;
        end
    end
end
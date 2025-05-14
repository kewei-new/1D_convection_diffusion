function error = evaluate_L2_norm_error(solution,X_points)

global nx hx
global weight

error = 0;
for i = 1:nx

    error = error + hx*weight'*(solution(:,i)-exact_fun(X_points(:,i))).^2;

end

error = sqrt(error);
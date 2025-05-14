function error = evaluate_absolute_norm_error(solution,X_points)

global nx
error = 0;
for i = 1:nx

    temp = max(abs(solution(:,i)-exact_fun(X_points(:,i))));
    error = max(error,temp);

end
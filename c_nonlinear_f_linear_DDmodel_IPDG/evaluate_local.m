function result = evaluate_local(x,uh,mp,derivate)
% evaluate the approximation in local segement
result = 0;

for m = 0:mp

    result = result + uh(m+1)*reference_basis(x,m,101,derivate);

end


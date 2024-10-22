function result = evaluate_ref(x,uh,mp,derivate)
% evaluate the approximation in local segement
% The input x needs to be radiated first to the reference cell [-0.5,0.5] by affine transformation
result = 0;

for m = 0:mp

    result = result + uh(m+1)*reference_basis(x,m,101,derivate);

end


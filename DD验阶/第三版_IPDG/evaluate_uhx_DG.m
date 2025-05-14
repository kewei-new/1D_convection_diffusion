function result = evaluate_uhx_DG(x,uh,mp)

result = 0;
for m = 0:mp

    result = result + uh(m+1)*reference_basis_x(x,m);

end

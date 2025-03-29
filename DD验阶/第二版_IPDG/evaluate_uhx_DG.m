function result = evaluate_uhx_DG(x,uh,mp)

result = 0;
ng = size(uh,2);
for n = 1:ng

    for m = 0:mp

    result = result + uh(m+1,n)*reference_basis_x(x,m);


    end
end
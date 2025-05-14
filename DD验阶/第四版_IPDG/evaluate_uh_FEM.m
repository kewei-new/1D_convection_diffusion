function result = evaluate_uh_FEM(x,uh,derivative,h)
% local

result = 0;
mp = size(uh,1);
mp =mp-1;
if mp == 1
    basis_type = 101;
else
    basis_type = 102;
end

for m = 1:mp+1

    result = result + uh(m)*reference_basis_FEM(x,derivative,m,basis_type)*h^(-derivative);

end
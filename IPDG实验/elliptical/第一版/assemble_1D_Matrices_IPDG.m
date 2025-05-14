function Matrices = assemble_1D_Matrices_IPDG(Gauss_coefficient,ng,mp)

Matrices = sparse(ng*(mp+1),ng*(mp+1));
Gp = Gauss_coefficient(:,1);
Gw = Gauss_coefficient(:,2);
for n = 1:ng


end
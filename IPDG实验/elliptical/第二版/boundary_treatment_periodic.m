function Matrices = boundary_treatment_periodic(Matrices)

global nx dimPk
% n = 1
Matrices(1:dimPk,(nx-1)*dimPk+1:nx*dimPk) = Matrices(dimPk+1:2*dimPk,1:dimPk);
% n = ng
Matrices((nx-1)*dimPk+1:nx*dimPk,1:dimPk) = Matrices(1:dimPk,dimPk+1:2*dimPk);
function [A,b]=boundary_treatment_FEM(boundary_nodes,A,b,Pb)

nbn = length(boundary_nodes);
for k = 1:nbn
    if boundary_nodes(1,k)==1
        % i：Finite element node index
        i = boundary_nodes(2,k);
        A(i,:) = 0;
        A(i,i) = 1;
        b(i) = boundary_fun(Pb(i));
    end
end

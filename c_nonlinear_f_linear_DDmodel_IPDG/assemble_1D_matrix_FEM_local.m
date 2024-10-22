function result = assemble_1D_matrix_FEM_local(M,P,T,Tb,mp)

ng = length(T);
result = sparse(ng*mp+1,ng*mp+1);

% Calculate the length of the interval
h_local = P(T(2,:)) - P(T(1,:));

for n = 1:ng
    for i = 1:mp+1
        for j = 1:mp+1

            result(Tb(i,n),Tb(j,n)) = result(Tb(i,n),Tb(j,n)) + M(i,j)/h_local(n);

        end
    end
end

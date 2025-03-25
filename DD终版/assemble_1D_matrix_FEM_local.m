function result = assemble_1D_matrix_FEM_local(Tb,h,basis_type)

ng = length(Tb);

if basis_type == 101
    m_FEM= 1;
    M = [1,-1;-1,1];
else
    m_FEM = 2;
    M = [7/3,-8/3,1/3;-8/3,16/3,-8/3;1/3,-8/3,7/3];
end
matrix_size = ng*m_FEM+1;
result = zeros(matrix_size,matrix_size);

for n = 1:ng
    for i = 1:m_FEM+1
        for j = 1:m_FEM+1

            result(Tb(i,n),Tb(j,n)) = result(Tb(i,n),Tb(j,n)) + M(i,j);

        end
    end
end

result = sparse(result/h);
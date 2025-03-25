function result = assemble_1D_matrix_with_variable(ng,mp,M2,L,H,R)
% 使用RK方法，每个时间步上都需要计算f(t+dt/2)，因此如果自变系数与时间有关就需要更新总纲矩阵
% 例如DD模型中(\mu*E*n)_x中E是与时间有关的量，因此每个时间半阶都需要先算出E，更新总纲矩阵
% update/generate matrix with variable coefficient
% M2: the matrice involing variable coefficient information
result = sparse(ng*(mp+1),ng*(mp+1));

for n = 2:ng-1
    result((n-1)*(mp+1)+1:n*(mp+1),(n-1)*(mp+1)+1:n*(mp+1)) =-M2(:,:,n)+H(:,:,n);
    result((n-1)*(mp+1)+1:n*(mp+1),n*(mp+1)+1:(n+1)*(mp+1)) =R(:,:,n+1);
    result((n-1)*(mp+1)+1:n*(mp+1),(n-2)*(mp+1)+1:(n-1)*(mp+1)) = L(:,:,n-1);
end

n=1;

result((n-1)*(mp+1)+1:n*(mp+1),(n-1)*(mp+1)+1:n*(mp+1)) =-M2(:,:,n)+H(:,:,n);
result((n-1)*(mp+1)+1:n*(mp+1),n*(mp+1)+1:(n+1)*(mp+1)) =R(:,:,n+1);
result((n-1)*(mp+1)+1:n*(mp+1),(ng-1)*(mp+1)+1:ng*(mp+1)) = L(:,:,ng);

n = ng;

result((n-1)*(mp+1)+1:n*(mp+1),(n-1)*(mp+1)+1:n*(mp+1)) =-M2(:,:,n)+H(:,:,n);
result((n-1)*(mp+1)+1:n*(mp+1),1:(mp+1)) =R(:,:,1);
result((n-1)*(mp+1)+1:n*(mp+1),(n-2)*(mp+1)+1:(n-1)*(mp+1)) = L(:,:,n);

end
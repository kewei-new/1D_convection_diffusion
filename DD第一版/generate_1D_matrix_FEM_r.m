function result = generate_1D_matrix_FEM_r(Gauss_coefficient,mp,trial_basis_derivative,test_basis_derivative,basis_type)

Gp = Gauss_coefficient(:,1);
Gw = Gauss_coefficient(:,2);
Gpn = length(Gw);

% num of basis fun
mpn = mp+1;
result = zeros(mpn,mpn);

for i = 1:mpn
    for j = 1:mpn
        for k = 1:Gpn
            
            result(i,j) = result(i,j) + Gw(k)*basis_reference_FEM(Gp(k),trial_basis_derivative,j,basis_type)...
                *basis_reference_FEM(Gp(k),test_basis_derivative,i,basis_type);

        end
    end
end

% 由于有限元在网格节点上的基函数是包含左右两部分的
% 因此在组装矩阵，在处理边界基函数内积的时候，需要按照左右两部分进行求和
% 其余基函数因为在边界处的函数值是0，因此只包含一部分就可以实现连续
% 不在这个函数里进行左右求和处理，是因为如果左右两部分是在不同区间内求和，如果区间是不均匀的就没法在参考单元上实现
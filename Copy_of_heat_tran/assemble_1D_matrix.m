function M = assemble_1D_matrix(Gauss_coefficient,inv_a,ng,mp,h,trial_drivative,test_derivative)

M = sparse(ng*(mp+1),ng*(mp+1));
GP = Gauss_coefficient(:,1);
Gw = Gauss_coefficient(:,2);

for n = 2:ng


    index_L = (n-2)*(mp+1);
    index = (n-1)*(mp+1);

    for alpha = 0:mp % test
        for beta = 0:mp % trial

            % element n
            index_1 = index + alpha+1;
            index_2 = index + beta+1;
            M(index_1,index_2) = inv_a(alpha+1,alpha+1)*(Gw'*(reference_basis(GP,beta,trial_drivative).*reference_basis(GP,alpha,test_derivative))*h^(1-trial_drivative-test_derivative)...
                               - reference_basis(0.5,beta,0)*reference_basis(0.5,alpha,0));

            % element n-1
            index_2 = index_L + beta+1;
            M(index_1,index_2) =inv_a(alpha+1,alpha+1)*(reference_basis(0.5,beta,0)*reference_basis(-0.5,alpha,0));


        end
    end
end

n = 1;
index_L = (ng-1)*(mp+1);
index = (n-1)*(mp+1);

for alpha = 0:mp % test
    for beta = 0:mp % trial

        % element n
        index_1 = index + alpha+1;
        index_2 = index + beta+1;
        M(index_1,index_2) = inv_a(alpha+1,alpha+1)*(Gw'*(reference_basis(GP,beta,trial_drivative).*reference_basis(GP,alpha,test_derivative))*h^(1-trial_drivative-test_derivative)...
                           - reference_basis(0.5,beta,0)*reference_basis(0.5,alpha,0));

        % element n-1
        index_2 = index_L + beta+1;
        M(index_1,index_2) = inv_a(alpha+1,alpha+1)*(reference_basis(0.5,beta,0)*reference_basis(-0.5,alpha,0));


    end
end

M = 2*M*h^(-1);


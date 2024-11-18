function M = assemble_1D_matrix_from_integral(Gauss_coefficient,inv_a,h,mp,trial_derivative,test_deribative)

M = zeros(mp+1,mp+1);
Gp = Gauss_coefficient(:,1);
Gw = Gauss_coefficient(:,2);

for alpha = 0:mp
    % test
    for beta = 0:mp
        % trial
        M(alpha+1,beta+1) = inv_a(alpha+1)*Gw'*(reference_basis(Gp,beta,trial_derivative).*reference_basis(Gp,alpha,test_deribative));
    end
end

M = (1/h)^(trial_derivative+test_deribative)*M;
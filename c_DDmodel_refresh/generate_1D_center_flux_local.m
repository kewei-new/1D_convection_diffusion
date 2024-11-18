function [H,R,L] = generate_1D_center_flux_local(inv_a,h,mp,trial_derivative,test_derivative,cha)
% when flux u = 1/2(u^+ + u^-)
H = zeros(mp+1,mp+1);
R = zeros(mp+1,mp+1);
L = zeros(mp+1,mp+1);

for alpha = 0:mp
    % test
    for beta = 0:mp
        % trial
        H(alpha+1,beta+1) = inv_a(alpha+1)*(1/2*reference_basis(0.5,beta,trial_derivative)*reference_basis(0.5,alpha,test_derivative)...
                      -1/2*reference_basis(-0.5,beta,trial_derivative)*reference_basis(-0.5,alpha,test_derivative) ...
                  +cha*(1/2*reference_basis(0.5,beta,test_derivative)*reference_basis(0.5,alpha,trial_derivative) ...
                      -1/2*reference_basis(-0.5,beta,test_derivative)*reference_basis(-0.5,alpha,trial_derivative)));
        R(alpha+1,beta+1) = inv_a(alpha+1)*(1/2*reference_basis(-0.5,beta,trial_derivative)*reference_basis(0.5,alpha,test_derivative) ...
                  -cha*1/2*reference_basis(-0.5,beta,test_derivative)*reference_basis(0.5,alpha,trial_derivative));
        L(alpha+1,beta+1) = inv_a(alpha+1)*((-1/2)*reference_basis(0.5,beta,trial_derivative)*reference_basis(-0.5,alpha,test_derivative) ...
                  +cha*1/2*reference_basis(0.5,beta,test_derivative)*reference_basis(-0.5,alpha,trial_derivative));

    end
end

H = (1/h)^(trial_derivative+test_derivative+1)*H;
R = (1/h)^(trial_derivative+test_derivative+1)*R;
L = (1/h)^(trial_derivative+test_derivative+1)*L;
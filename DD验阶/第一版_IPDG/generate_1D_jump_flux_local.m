function [H,R,L] = generate_1D_jump_flux_local(inv_a,h,mp)

H = zeros(mp+1,mp+1);
R = zeros(mp+1,mp+1);
L = zeros(mp+1,mp+1);

for alpha = 0:mp
    % test
    for beta = 0:mp
        % trial
        H(alpha+1,beta+1) = inv_a(alpha+1)*(reference_basis(0.5,beta)*reference_basis(0.5,alpha)...
                      +reference_basis(-0.5,beta)*reference_basis(-0.5,alpha));
        R(alpha+1,beta+1) = inv_a(alpha+1)*(-reference_basis(-0.5,beta)*reference_basis(0.5,alpha));
        L(alpha+1,beta+1) = inv_a(alpha+1)*(-reference_basis(0.5,beta)*reference_basis(-0.5,alpha));

    end
end

H = (1/h)*H;
R = (1/h)*R;
L = (1/h)*L;
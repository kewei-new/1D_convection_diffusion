function [H,R,L] = generate_1D_center_flux_local(inv_a,h,mp,cha)
% when flux u = 1/2(u^+ + u^-)
H = zeros(mp+1,mp+1);
R = zeros(mp+1,mp+1);
L = zeros(mp+1,mp+1);

for alpha = 0:mp
    % test
    for beta = 0:mp
        % trial
        H(alpha+1,beta+1) = inv_a(alpha+1)*(1/2*reference_basis_x(0.5,beta)*reference_basis(0.5,alpha)...
                      -1/2*reference_basis_x(-0.5,beta)*reference_basis(-0.5,alpha) ...
                  +cha*(1/2*reference_basis(0.5,beta)*reference_basis_x(0.5,alpha) ...
                      -1/2*reference_basis(-0.5,beta)*reference_basis_x(-0.5,alpha)));
        R(alpha+1,beta+1) = inv_a(alpha+1)*(1/2*reference_basis_x(-0.5,beta)*reference_basis(0.5,alpha) ...
                  -cha*1/2*reference_basis(-0.5,beta)*reference_basis_x(0.5,alpha));
        L(alpha+1,beta+1) = inv_a(alpha+1)*((-1/2)*reference_basis_x(0.5,beta)*reference_basis(-0.5,alpha) ...
                  +cha*1/2*reference_basis(0.5,beta)*reference_basis_x(-0.5,alpha));

    end
end

H = (1/h)^(2)*H;
R = (1/h)^(2)*R;
L = (1/h)^(2)*L;
function [H,R,L] = generate_1D_penalty_flux(inv_a,mp,h)

H = zeros(mp+1,mp+1);
R = zeros(mp+1,mp+1);
L = zeros(mp+1,mp+1);

for m1 = 0:mp
    % test

    for m2 = 0:mp
        % trial

        H(m1+1,m2+1) = reference_basis(0.5,m2,0)*reference_basis(0.5,m1,0)...
                      +reference_basis(-0.5,m2,0)*reference_basis(-0.5,m1,0);
        R(m1+1,m2+1) = -reference_basis(-0.5,m2,0)*reference_basis(0.5,m1,0);
        L(m1+1,m2+1) = -reference_basis(0.5,m2,0)*reference_basis(-0.5,m1,0);
    end

end

H = inv_a*H*h^(-1);
R = inv_a*R*h^(-1);
L = inv_a*L*h^(-1);

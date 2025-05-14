function generate_1D_constant_Matrices_IPDG(Gauss_coefficient,ng,h,mp)

Gp = Gauss_coefficient(:,1);
Gw = Gauss_coefficient(:,2);

M1 = zeros(mp+1,mp+1);
M2L = zeros(mp+1,mp+1,2);
M2R = zeros(mp+1,mp+1,2);
M3L = zeros(mp+1,mp+1,2);
M3R = zeros(mp+1,mp+1,2);

% (ux,vx)
for d1 = 0:mp
    for d2 = 0:mp

        M1(d1+1,d2+1) = Gw'*1/h*reference_basis_x(Gp,d2)*reference_basis_x(Gp,d1);

    end
end

% ({ux},v)e = 1/2*(ux(0.5) + ux(-0.5))v(0.5) - 1/2*(uxL(0.5) + ux(-0.5))v(-0.5)
% 0.5*(ux(0.5)v(0.5) + uxR(-0.5)v(0.5) - uxL(0.5)v(-0.5) - ux(-0.5)v(-0.5))
for d1 = 0:mp
    for d2 = 0:mp

        M2L(d1+1,d2+1,1) =-1/h*reference_basis_x(-0.5,d2)*reference_basis(-0.5,d1);
        M2L(d1+1,d2+1,2) =-1/h*reference_basis_x(0.5,d2)*reference_basis(-0.5,d1);
        M2R(d1+1,d2+1,1) = 1/h*reference_basis_x(0.5,d2)*reference_basis(0.5,d1);
        M2R(d1+1,d2+1,2) = 1/h*reference_basis_x(-0.5,d2)*reference_basis(0.5,d1);

    end
end

% 0.5*sign(du/dn)([u],vx)e = 1/2*(u(0.5-)-u(0.5+))vx(0.5) + 1/2*(u(-0.5-)-u(-0.5+))vx(-0.5)

% -([u],v)e =-(u(0.5-)-u(0.5+))v(0.5-) + (u(-0.5-)-u(-0.5+))v(-0.5+)
% uR(-0.5)v(0.5) - u(0.5)v(0.5) + uL(0.5)v(-0.5) - u(-0.5)v(-0.5)
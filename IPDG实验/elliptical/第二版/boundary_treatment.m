function [Matrices,Vector] = boundary_treatment(Matrices,Vector)
% Dirichlet B.C
global dimPk nx hx
global phixGL phixGR phiGL phiGR
global beta alpha

% (gD,(beta*dv/dn+alpha/h*v))

% left
for d = 1:dimPk

    Vector(d,1) = Vector(d,1) + hx*Dirichlet_fun(0)*(-beta*phixGL(:,d)-alpha/hx*phiGR(:,d));

end
% right
for d = 1:dimPk

    Vector((nx-1)*dimPk+d,1) = Vector((nx-1)*dimPk+d,1) + hx*Dirichlet_fun(2*pi)*(beta*phixGR(:,d)-alpha/hx*phiGL(:,d));

end

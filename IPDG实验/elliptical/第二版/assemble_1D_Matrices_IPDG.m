function Matrices = assemble_1D_Matrices_IPDG

global hx nx dimPk alpha beta
global phixG phiGL phiGR phixGL phixGR
global weight

M1 = zeros(dimPk,dimPk);M2LL = zeros(dimPk,dimPk);M2RR = zeros(dimPk,dimPk);M2RL = zeros(dimPk,dimPk);M2LR = zeros(dimPk,dimPk);
M3LL = zeros(dimPk,dimPk);M3RR = zeros(dimPk,dimPk);M3RL = zeros(dimPk,dimPk);M3LR = zeros(dimPk,dimPk);
% -(ux,vx) + ({ux},v)e + 0.5*beta*sign(du/dn)*([u],vx)e - alpha*([u],v)e

% (ux,vx)
for d1 = 1:dimPk %test
    for d2 = 1:dimPk %trial

        M1(d1,d2) = hx*weight'*(phixG(:,d2).*phixG(:,d1));

    end
end

% ({ux},v)e -> 0.5(ux,v)L & 0.5(ux,v)R
for d1 = 1:dimPk %test
    for d2 = 1:dimPk %trial

        M2LL(d1,d2) = phixGL(:,d2).*phiGL(:,d1);
        M2RR(d1,d2) = phixGR(:,d2).*phiGR(:,d1);
        M2RL(d1,d2) = phixGR(:,d2).*phiGL(:,d1);
        M2LR(d1,d2) = phixGL(:,d2).*phiGR(:,d1);

    end
end
M2LL = 0.5*M2LL;
M2RR = 0.5*M2RR;
M2RL = 0.5*M2RL;
M2LR = 0.5*M2LR;

% ([u],vx)e -> (u,vx)L & (u,vx)R
% (u,vx)L = M2L'; (u,vx)R = M2R';

% alpha/hx*([u],v)e -> alpha/hx*(u,v)L & alpha/hx*(u,v)R
for d1 = 1:dimPk %test
    for d2 = 1:dimPk %trial

        M3LL(d1,d2) = phiGL(:,d2).*phiGL(:,d1);
        M3RR(d1,d2) = phiGR(:,d2).*phiGR(:,d1);
        M3RL(d1,d2) = phiGR(:,d2).*phiGL(:,d1);
        M3LR(d1,d2) = phiGL(:,d2).*phiGR(:,d1);

    end
end
M3LL = alpha/hx*M3LL;
M3RR = alpha/hx*M3RR;
M3RL = alpha/hx*M3RL;
M3LR = alpha/hx*M3LR;

% assemble general Matrices
Matrices = sparse(nx*dimPk,nx*dimPk);
index = @(i)(i-1)*dimPk;
for i = 1:nx

    Matrices(index(i)+1:index(i+1),index(i)+1:index(i+1)) = -M1 + (M2RR-M2LL) + beta*(M2RR'-M2LL') - (M3RR+M3LL);
    if i > 1
        Matrices(index(i)+1:index(i+1),index(i-1)+1:index(i)) =-M2RL + beta*M2LR' + M3RL;
    end
    if  i < nx
        Matrices(index(i)+1:index(i+1),index(i+1)+1:index(i+2)) = M2LR - beta*M2RL' + M3LR;
    end

end

% 如果是周期边界，需要在这之后再加上边界处理
% Matrices = boundary_treatment_periodic(Matrices);

function get_basis

global hx dimPk
global lambda
global phiG phixG phiGL phiGR phixGL phixGR phimass

NumGLP = length(lambda);
phiG = zeros(NumGLP,4);phiGL = zeros(1,4);phiGR = zeros(1,4);
phixG = zeros(NumGLP,4);phixGL = zeros(1,4);phixGR = zeros(1,4);
phimass = zeros(4,1);

for i1 = 1:NumGLP
    % 1
    phiG(i1,1) = 1;
    phiGL(1,1) = 1;
    phiGR(1,1) = 1;
    phixG(i1,1) = 0;
    phixGL(1,1) = 0;
    phixGR(1,1) = 0;
    phimass(1) = 1;

    % x
    phiG(i1,2) = lambda(i1);
    phiGL(1,2) = -0.5;
    phiGR(1,2) = 0.5;
    phixG(i1,2) = 1/hx;
    phixGL(1,2) = 1/hx;
    phixGR(1,2) = 1/hx;
    phimass(2) = 1/12.0;
    
    % x.^2 - 1./12;
    phiG(i1,3) = lambda(i1)^2 - 1/12.0;
    phiGL(1,3) = 1/6.0;
    phiGR(1,3) = 1/6.0;
    phixG(i1,3) = 2*lambda(i1)/hx;
    phixGL(1,3) =-1/hx;
    phixGR(1,3) = 1/hx;
    phimass(3) = 1/180.0;

    % x.^3 - 0.15*x;
    phiG(i1,4) = lambda(i1)^3 - 0.15*lambda(i1);
    phiGL(1,4) =-0.0500;
    phiGR(1,4) = 0.0500;
    phixG(i1,4) = (3*lambda(i1)^2 - 0.15)/hx;
    phixGL(1,4) = 0.6/hx;
    phixGR(1,4) = 0.6/hx;
    phimass(4) = 1/2800.0;
    
end

phiG = phiG(:,1:dimPk);
phixG = phixG(:,1:dimPk);
phiGL = phiGL(1,1:dimPk);
phiGR = phiGR(1,1:dimPk);
phixGL = phixGL(1,1:dimPk);
phixGR = phixGR(1,1:dimPk);
phimass = phimass(1:dimPk);
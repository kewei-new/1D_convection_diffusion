function gauss = generate_1D_gaussian_quadrature_r(mp)
if mp == 0
    gauss(1,1) = -0.5;
    gauss(2,1) = 0.5;
    gauss(1,2) = 0.5;
    gauss(2,2) = 0.5;
end

if mp == 1
    % the points of 4th order Gauss-Lobatto quadrature
    gauss(1,1) = -0.5;
    gauss(2,1) = 0.5;
    gauss(3,1) = 0.0;
    % coefficients of 4th order Gauss-Lobatto quadrature
    gauss(1,2) = 1.0/6.0;
    gauss(2,2) = gauss(1,2);
    gauss(3,2) = 2.0/3.0;
end

if mp == 2
    % the points of 6th order Gauss-Lobatto quadrature
    gauss(1,1) = -0.5;
    gauss(2,1) = 0.5;
    gauss(3,1) = -sqrt(5)/10.0;
    gauss(4,1) = sqrt(5)/10.0;
    % coefficients of 6th order Gauss-Lobatto quadrature
    gauss(1,2) = 1.0/12.0;
    gauss(2,2) = gauss(1,2);
    gauss(3,2) = 5.0/12.0;
    gauss(4,2) = gauss(3,2);
end

if mp == 3
    % the points of 8th order Gauss-Lobatto quadrature
    gauss(1,1) = -0.5;
    gauss(2,1) = 0.5;
    gauss(3,1) = -sqrt(21)/14.0;
    gauss(4,1) = sqrt(21)/14.0;
    gauss(5,1) = 0.0;
    % coefficients of 8th order Gauss-Lobatto quadrature
    gauss(1,2) = 1.0/20.0;
    gauss(2,2) = gauss(1,2);
    gauss(3,2) = 49.0/180.0;
    gauss(4,2) = gauss(3,2);
    gauss(5,2) = 64.0/180.0;
end

if mp == 4
    % the points of 10th order Gauss-Lobatto quadrature
    gauss(1,1) = -0.5;
    gauss(2,1) = 0.5;
    gauss(3,1) = -sqrt(147 + 42*sqrt(7))/42.0;
    gauss(4,1) = sqrt(147 + 42*sqrt(7))/42.0;
    gauss(5,1) = -sqrt(147 - 42*sqrt(7))/42.0;
    gauss(6,1) = sqrt(147 - 42*sqrt(7))/42.0;
    % coefficients of 10th order Gauss-Lobatto quadrature
    gauss(1,2) = 1.0/30.0;
    gauss(2,2) = gauss(1,2);
    gauss(3,2) = (-7 + 5*sqrt(7))*sqrt(7)*(7 + sqrt(7))/840.0;
    gauss(4,2) = gauss(3,2);
    gauss(5,2) = (7 + 5*sqrt(7))*sqrt(7)/(7 + sqrt(7))/20.0;
    gauss(6,2) = gauss(5,2);
end

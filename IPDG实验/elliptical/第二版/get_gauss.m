function [lambda,weight] = get_gauss(mp)

if mp == 1

    lambda = zeros(3,1);
    weight = zeros(3,1);
    % the points of 4th order Gauss-Lobatto quadrature
    lambda(1) = -0.5;
    lambda(2) = 0.5;
    lambda(3) = 0;
    % coefficients of 4th order Gauss-Lobatto quadrature
    weight(1) = 1/6.0;
    weight(2) = 1/6.0;
    weight(3) = 2/3.0;
elseif mp == 2

    lambda = zeros(4,1);
    weight = zeros(4,1);
    % the points of 6th order Gauss-Lobatto quadrature
    lambda(1) = -0.5;
    lambda(2) = 0.5;
    lambda(3) = -sqrt(5)/10.0;
    lambda(4) = sqrt(5)/10.0;
    % coefficients of 6th order Gauss-Lobatto quadrature
    weight(1) = 1.0/12.0;
    weight(2) = weight(1);
    weight(3) = 5.0/12.0;
    weight(4) = weight(3);

elseif mp == 3

    lambda = zeros(5,1);
    weight = zeros(5,1);    
    % the points of 8th order Gauss-Lobatto quadrature
    lambda(1) = -0.5;
    lambda(2) = 0.5;
    lambda(3) = -sqrt(21)/14.0;
    lambda(4) = sqrt(21)/14.0;
    lambda(5) = 0.0;
    % coefficients of 8th order Gauss-Lobatto quadrature
    weight(1) = 1.0/20.0;
    weight(2) = weight(1);
    weight(3) = 49.0/180.0;
    weight(4) = weight(3);
    weight(5) = 64.0/180.0;
else

    lambda = zeros(6,1);
    weight = zeros(6,1);     
    % the points of 10th order Gauss-Lobatto quadrature
    lambda(1) = -0.5;
    lambda(2) = 0.5;
    lambda(3) = -sqrt(147 + 42*sqrt(7))/42.0;
    lambda(4) = sqrt(147 + 42*sqrt(7))/42.0;
    lambda(5) = -sqrt(147 - 42*sqrt(7))/42.0;
    lambda(6) = sqrt(147 - 42*sqrt(7))/42.0;
    % coefficients of 10th order Gauss-Lobatto quadrature
    weight(1) = 1.0/30.0;
    weight(2) = weight(1);
    weight(3) = (-7 + 5*sqrt(7))*sqrt(7)*(7 + sqrt(7))/840.0;
    weight(4) = weight(3);
    weight(5) = (7 + 5*sqrt(7))*sqrt(7)/(7 + sqrt(7))/20.0;
    weight(6) = weight(5);

end

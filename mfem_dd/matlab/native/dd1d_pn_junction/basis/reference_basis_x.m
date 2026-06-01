function res = reference_basis_x(x,mp)
% x:coordinate
% mp:hightest order of the polynomial space
% basis_derivate:the derivate order of the basis funciton is 1


if mp == 0
    res = 0.*x;

elseif mp == 1
    res = 0.*x + 1;

elseif mp == 2
    res = 2.*x;

elseif mp == 3
    res = 3*x.^2 - 0.15;

elseif mp == 4
    res = 4*x.^3 - 3/7*x;

end
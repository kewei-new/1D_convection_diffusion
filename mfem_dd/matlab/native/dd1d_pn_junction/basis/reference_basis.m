function res = reference_basis(x,mp)
% x:coordinate
% mp:hightest order of the polynomial space
% basis_derivate:the derivate order of the basis funciton is 0

if mp == 0
    res = 0.*x+1;

elseif mp == 1
    res = 1.*x;

elseif mp == 2
    res = x.^2 - 1./12;

elseif mp == 3
    res = x.^3 - 0.15*x;

elseif mp == 4
    res = (x.^2-3/14).*x.*x+3.0/560;

end





end
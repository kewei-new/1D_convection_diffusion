function [hx,Xm] = generate_1D_uniform_grid(xa,xb)

global nx
hx = (xb-xa)/nx;
Xm = zeros(nx+1,1);
for i = 1:nx

    Xm(i) = xa + hx*(i-1/2);

end
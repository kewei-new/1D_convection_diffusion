function [P,T] = generate_grid_1D(left,right,ng)

P = linspace(left,right,ng+1);
for n = 1:ng
    T(1,n) = n;
    T(2,n) = n+1;
end
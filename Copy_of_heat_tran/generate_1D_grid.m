function [P,T] = generate_1D_grid(left,right,ng)

P = linspace(left,right,ng+1);
for n = 1:ng
    T(1,n) = n;
    T(2,n) = n+1;
end
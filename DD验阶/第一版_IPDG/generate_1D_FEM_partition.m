function [Pb,Tb] = generate_1D_FEM_partition(P,T)

mp = 2;

ng = length(T);
left = P(1);
right = P(ng+1);

Pb = linspace(left,right,1+mp*ng); % ng+1+(mp-1)*ng
Tb = zeros(mp+1,ng);% 2+(mp-1)

for n = 1:ng

    for ns = 1:mp+1

         Tb(ns,n) = ns + (n-1)*mp;

    end
end

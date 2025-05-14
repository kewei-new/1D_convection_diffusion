function [h,mid_points] = generate_1D_uneven_grid(left,right,ng)

temp = zeros(1,ng+1);
mid_points = zeros(1,ng);
h = zeros(1,ng);

n1 = ng/4;
n2 = ng/2;
n3 = ng/4;
dx1 = (right-left)/(3*n1);
dx2 = (right-left)/(3*n2);
dx3 = (right-left)/(3*n3);

temp(1) = left;
for n = 1:ng
    
    if n<=n1
        temp(n+1) = temp(1) + dx1*n;
        h(n) = dx1;

    elseif n<=n2+n1 && n>n1
        temp(n+1) = temp(n1+1) + dx2*(n-n1);
        h(n) = dx2;
    else
        temp(n+1) = temp(n2+n1+1) + dx3*(n-n1-n2);
        h(n) = dx3;
    end

    mid_points(n) = (temp(n)+temp(n+1))/2;
end
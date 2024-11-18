function b = generate_1D_vector(t,Gauss_coefficient,inv_a,mid_points,mp,ng,h)

b = zeros(mp+1,ng);
Gp = Gauss_coefficient(:,1);
Gw = Gauss_coefficient(:,2);
Gpn = length(Gw);

local_Gp = Gp*h + mid_points;
local_Gp = reshape(local_Gp,size(local_Gp,1)*size(local_Gp,2),1);
f_value = f_fun1(local_Gp,t);
f_value = repmat(f_value,1,mp+1);

v_value = zeros(Gpn,mp+1);
for m = 0:mp
    v_value(:,m+1) = inv_a(m+1,m+1)*reference_basis(Gp,m,0);
end

v_value = repmat(v_value,ng,1);
% size(fv_value) == (Gpn*ng,mp+1)
% 第m列表示f(Gp)*vm(Gp)
fv_value = f_value.*v_value; 

for n = 1:ng
    for m = 0:mp

        b(m+1,n) = Gw'*fv_value((n-1)*Gpn+1:n*Gpn,m+1);

    end
end

function Vector = generate_1D_vector

global nx hx Xm dimPk
global phiG
global lambda weight

Vector = zeros(nx*dimPk,1);
index = @(i)(i-1)*dimPk;
for i = 1:nx
    f_value = f_fun(lambda(:)*hx+Xm(i));
    for d = 1:dimPk  
        Vector(index(i)+d,1) = hx*weight'*(f_value.*phiG(:,d));
    end
end
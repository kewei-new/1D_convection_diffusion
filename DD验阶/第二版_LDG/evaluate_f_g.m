function result = evaluate_f_g(Gauss_coefficient,inv_mass,phixh,uh1,mp,h,mid_points,emf,t)

ng = size(uh1,2);
f = zeros(mp+1,ng);
g = zeros(mp+1,ng);
hu = zeros(mp+1,ng);
Gp = Gauss_coefficient(:,1);
Gw = Gauss_coefficient(:,2);


u1 = ones(length(Gp),ng+2);
phix = ones(length(Gp),ng+2);
for k = 2:ng+1

   u1(:,k) = evaluate_uh_DG(Gp,uh1(:,k-1),mp);
   phix(:,k) = evaluate_uh_DG(Gp,phixh(:,k),mp);
   
end

u1(2,1) = u1(2,ng+1);
phix(2,1) = phix(2,ng+1);

u1(1,ng+2) = u1(1,2);
phix(1,ng+2) = phix(1,2);

% u1 = n
n = u1;

% 已检查
% evaluate -(f1,vx) + (f1,v)_e
% f1 = -phix.*n

f1 = n.*(-phix);

for k = 2:ng+1
    for m = 0:mp
    
        f(m+1,k-1) = inv_mass(m+1)*(-Gw'*(f1(:,k).*reference_basis_x(Gp,m))...
                  +1/2*((f1(2,k)+f1(1,k+1))*reference_basis(0.5,m)-(f1(2,k-1)+f1(1,k))*reference_basis(-0.5,m))...
                  -1/2*emf*((u1(1,k+1)-u1(2,k))*reference_basis(0.5,m)-(u1(1,k)-u1(2,k-1))*reference_basis(-0.5,m)))*h(k-1)^(-1);


    end
end

% evaluate (g1xx,v)
% g1 = n
qh = zeros(mp+1,ng);
q1 = zeros(length(Gp),ng+1);

g1 = n;

for k = 2:ng+1
    for m = 0:mp
    
        qh(m+1,k-1) =inv_mass(m+1)*(-Gw'*(g1(:,k).*reference_basis_x(Gp,m))...
              +g1(2,k)*reference_basis(0.5,m)-g1(2,k-1)*reference_basis(-0.5,m))*h(k-1)^(-1);
    
    end
    q1(:,k-1) = evaluate_uh_DG(Gp,qh(:,k-1),mp);
end
q1(1,ng+1) = q1(1,1);

for k =1:ng
    for m = 0:mp
    
        g(m+1,k) =inv_mass(m+1)*(-Gw'*(q1(:,k).*reference_basis_x(Gp,m))...
              -q1(1,k)*reference_basis(-0.5,m)+q1(1,k+1)*reference_basis(0.5,m))*h(k)^(-1);
    
    end
end


for k = 2:ng+1
    for m = 0:mp
    
        hu(m+1,k-1) = inv_mass(m+1)*(Gw'*(f1_fun(Gp*h(k-1) + mid_points(k-1),t).*reference_basis(Gp,m)));
    
    end
end


result = g + f + hu;

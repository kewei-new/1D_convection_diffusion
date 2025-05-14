function result = evaluate_f_g(Gauss_coefficient,inv_mass,Tb,phih,uh1,mp,h,mid_points,emf,C_div_h,t)

% 控制
method = 1;

ng = size(uh1,2) - 2;
f = zeros(mp+1,ng+2);
g = zeros(mp+1,ng+2);
Gp = Gauss_coefficient(:,1);
Gw = Gauss_coefficient(:,2);


u1 = ones(length(Gp),ng+2);
u1x = ones(length(Gp),ng+2);
E = ones(length(Gp),ng+2);

phih = phih(Tb);
for k = 2:ng+1

   u1(:,k) = evaluate_uh_DG(Gp,uh1(:,k),mp);
   u1x(:,k) = evaluate_uhx_DG(Gp,uh1(:,k),mp)*h^(-1);
   E(:,k) = evaluate_uh_FEM(Gp,phih(:,k),1,h);
   
end

u1(2,1) = u1(2,ng+1);
u1x(2,1) = u1x(2,ng+1);
E(2,1) = E(2,ng+1);

u1(1,ng+2) = u1(1,2);
u1x(1,ng+2) = u1x(1,2);
E(1,ng+2) = E(1,2);

% 已检查
% evaluate -(f1,vx) + (f1,v)_e
% f1 = -phix.*n

f1 = u1.*E;

for k = 2:ng+1
    for m = 0:mp
    
        f(m+1,k,1) = inv_mass(m+1)*(-Gw'*(f1(:,k).*reference_basis_x(Gp,m))...
                  +1/2*((f1(2,k)+f1(1,k+1))*reference_basis(0.5,m)-(f1(2,k-1)+f1(1,k))*reference_basis(-0.5,m))...
                  -1/2*emf*((u1(1,k+1)-u1(2,k))*reference_basis(0.5,m)-(u1(1,k)-u1(2,k-1))*reference_basis(-0.5,m)))*h^(-1);

            % f(m+1,k) = inv_mass(m+1)*(-Gw'*(f1(:,k).*reference_basis_x(Gp,m))...
            %   +1/2*((f1(2,k)+f1(1,k+1))*reference_basis(0.5,m)-(f1(2,k-1)+f1(1,k))*reference_basis(-0.5,m))...
            %   )*h^(-1);

    end
end

% evaluate -(g1x,vx) + (g1x,v)e
% g1 = DE.*n

g1x = u1x;
g1 = u1;

g11 =zeros(mp+1,ng+2);

for k = 2:ng+1
    for m = 0:mp

        % g11(m+1,k) =  inv_mass(m+1)*(-1/2*g1x(2,k-1)*reference_basis(-0.5,m)...
        %       +1/2*g1(2,k-1)*reference_basis_x(-0.5,m)*h^(-1)*method)*h^(-1);
        % 
        g(m+1,k) =inv_mass(m+1)*(-Gw'*(g1x(:,k).*reference_basis_x(Gp,m))...
              +1/2*((g1x(2,k)+g1x(1,k+1))*reference_basis(0.5,m)-(g1x(2,k-1)+g1x(1,k))*reference_basis(-0.5,m))...
              +1/2*((g1(2,k)-g1(1,k+1))*reference_basis_x(0.5,m)+(g1(2,k-1)-g1(1,k))*reference_basis_x(-0.5,m))*h^(-1)*method...
              -C_div_h*((g1(2,k)-g1(1,k+1))*reference_basis(0.5,m)+(g1(1,k)-g1(2,k-1))*reference_basis(-0.5,m)))*h^(-1);
    
    end
end



hu = zeros(mp+1,ng+2);
for k = 2:ng+1
    for m = 0:mp
    
        hu(m+1,k) = inv_mass(m+1)*(Gw'*(f1_fun(Gp*h + mid_points(k-1),t).*reference_basis(Gp,m)));
    
    end
end


result = g - f + hu;

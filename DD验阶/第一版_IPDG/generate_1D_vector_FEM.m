function result=generate_1D_vector_FEM(Gauss_coefficient,t,uh,Tb,h,mid_points,mp,ng,basis_type)

if basis_type == 101
    m_FEM=1;
else
    m_FEM=2;
end

Gp = Gauss_coefficient(:,1);
Gw = Gauss_coefficient(:,2);
result = zeros(ng*m_FEM+1,1);

% nd
Gp_local = Gp*h + mid_points;
nd_value = nd(Gp_local,t);

% n
temp_value = zeros(length(Gp),mp+1);
for m = 0:mp  
    temp_value(:,m+1) = reference_basis(Gp,m); 
end
n_value = temp_value*uh;

% n-nd
f_value = -(n_value-nd_value);

% fv

for n = 1:ng
    
    for m = 1:m_FEM+1

        result(Tb(m,n),1) = result(Tb(m,n),1)+Gw'*(f_value(:,n).*reference_basis_FEM(Gp,0,m,basis_type))*h;
        
    end
end


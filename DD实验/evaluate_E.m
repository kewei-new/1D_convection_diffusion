function result = evaluate_E(phi,Tb,basis_type,ng,h)

result = zeros(2,ng);
phi = phi(Tb);

if basis_type == 101
    m_FEM = 1;
else
    m_FEM = 2;
end

for n = 1:ng

    for m = 1:m_FEM + 1
        
        result(1,n) = result(1,n) - phi(m,n)*reference_basis_FEM(-0.5,1,m,basis_type)*h^(-1);
        result(2,n) = result(2,n) - phi(m,n)*reference_basis_FEM( 0.5,1,m,basis_type)*h^(-1);

    end

end

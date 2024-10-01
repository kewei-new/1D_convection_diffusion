function temp(legendre_basis)
a = zeros(5,5);

for i = 0:4

    for j = i:4

        f = @(x)(legendre_basis(x,i,101,0).*legendre_basis(x,j,101,0));
        a(i+1,j+1) = integral(f,-0.5,0.5);

    end

end

a
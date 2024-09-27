function res = reference_basis(x,mp,basis_type,basis_derivative)
% x:coordinate
% mp:hightest order of the polynomial space
% basis_type:the type of the basis funciton
% basis_derivate:the derivate order of the basis funciton

if basis_type == 101

    if basis_derivative == 0

        if mp == 0
            res = 0.*x+1;   

        elseif mp == 1        
            res = 1.*x;  

        elseif mp == 2       
            res = x.^2 - 1./12;
        
        elseif mp == 3
            res = x.^3 - 0.15*x;
        
        elseif mp == 4
            res = (x.^2-3/14).*x.*x+3.0/560;
        
        end

    elseif basis_derivative == 1
        
        if mp == 0
            res = 0.*x;   

        elseif mp == 1        
            res = 0.*x + 1;  

        elseif mp == 2       
            res = 2.*x;
        
        elseif mp == 3
            res = 3*x.^2 - 0.15;
        
        elseif mp == 4
            res = 4*x.^3 - 3/7*x;
        
        end

    end

end




end
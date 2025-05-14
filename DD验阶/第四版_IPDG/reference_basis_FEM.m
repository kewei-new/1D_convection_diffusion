function result = reference_basis_FEM(x,basis_derivative,basis_index,basis_type)
% in [-0.5,0.5]

if basis_type == 101

    if basis_derivative == 0

       if basis_index == 1

           result = 0.5 - x;
    
       elseif basis_index == 2

           result = 0.5 + x;

       else
           error = "wrong index!"
       end

    elseif basis_derivative == 1
        
       if basis_index == 1

           result = -1 + 0.*x;

       elseif basis_index == 2
            
           result = 1 + 0.*x;
           
       else
           error = "wrong index!"
       end

    end

elseif basis_type == 102

    if basis_derivative == 0

        if basis_index == 1

            result = 2*x.^2 - x;

        elseif basis_index == 2

            result = -4*x.^2 + 1;

        elseif basis_index == 3

            result = 2*x.^2 + x;

        else
           error = "wrong index!"
        end

    elseif basis_derivative == 1
        
        if basis_index == 1

            result = 4*x - 1;

        elseif basis_index == 2

            result = -8*x;

        elseif basis_index == 3

            result = 4*x + 1;

        else
           error = "wrong index!"
        end
     
    elseif basis_derivative == 2
        
        if basis_index == 1

            result = 4 + 0.*x;

        elseif basis_index == 2

            result = -8 + 0.*x;

        elseif basis_index == 3

            result = 4 + 0.*x;

        else
           error = "wrong index!"
        end

    end

elseif basis_type == 103

    if basis_derivative == 0

        if basis_index == 1

            result = -4.5*x.^3 + 2.25*x.^2 + 0.125*x - 0.0625;

        elseif basis_index == 2

            result = 13.5*x.^3 - 2.25*x.^2 - 3.375*x + 0.5625;

        elseif basis_index == 3

            result = -13.5*x.^3 - 2.25*x.^2 - 3.375*x + 0.5625;

        elseif basis_index == 4

            result = 4.5*x.^3 + 2.25*x.^2 - 0.125*x - 0.0625;

        else
           error = "wrong index!"
        end

    elseif basis_derivative == 1
        
        if basis_index == 1

            result = -13.5*x.^2 + 4.5*x + 0.125;

        elseif basis_index == 2

            result = 40.5*x.^2 - 4.5*x -3.375;

        elseif basis_index == 3

            result = -40.5*x.^2 - 4.5*x -3.375;

        elseif basis_index == 4

            result = 13.5*x.^2 + 4.5*x - 0.125;

        else
           error = "wrong index!"
        end
     
    elseif basis_derivative == 2
        
        if basis_index == 1

        elseif basis_index == 2

        elseif basis_index == 3

        elseif basis_index == 4

        else
           error = "wrong index!"
        end

    elseif basis_derivative == 3
        
        if basis_index == 1

        elseif basis_index == 2

        elseif basis_index == 3

        elseif basis_index == 4

        else
           error = "wrong index!"
        end        

    end

end




end
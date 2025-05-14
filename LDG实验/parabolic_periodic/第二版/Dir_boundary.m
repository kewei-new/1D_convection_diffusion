function result = Dir_boundary(x,t)

if abs(x)<1e-10

    result = exp(-t)*sin(-t);
    % result = sin(-t);

elseif abs(x-2*pi)<1e-10
    
    result = exp(-t)*sin(-t);
    % result = sin(-t);

end


end
function result = LDG_bound(x,t)

if abs(x)<1e-10
    result = sin(t);

elseif abs(x-2*pi)<1e-10
    result = sin(t);

end

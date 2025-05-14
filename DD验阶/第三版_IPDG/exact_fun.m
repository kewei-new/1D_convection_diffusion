function result = exact_fun(x,t,i)

if i==1
    result = sin(x)*cos(t);
elseif i==2

    result = sin(t)*cos(x);
elseif i==3

    result = sin(t)*sin(x);
end
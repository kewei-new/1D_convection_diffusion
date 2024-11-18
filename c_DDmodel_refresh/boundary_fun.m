function result = boundary_fun(x)

global xk T_0 xe ni vbias

if x==0
    C = xk*T_0/xe;
    result = C*log(nd(x)/ni);

elseif x==0.6

    C = xk*T_0/xe;
    result = C*log(nd(x)/ni) + vbias;

end

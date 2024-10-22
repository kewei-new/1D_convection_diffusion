function result = boundary_fun(x)

load('boundary.mat');

if x==0
    C = boundary.xk*boundary.T_0/boundary.xe;
    result = C*log(nd(x)/boundary.ni);

elseif x==0.6

    C = boundary.xk*boundary.T_0/boundary.xe;
    result = C*log(nd(x)/boundary.ni) + boundary.vbias;

end

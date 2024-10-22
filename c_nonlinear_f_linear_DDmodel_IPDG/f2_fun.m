function result = f2_fun(x,x_local,uh_local)
% e/{\epsilon}*(n - nd)

load('params.mat');
C = params.xe/params.epsilon;
mp = length(uh_local)-1;

result = -C*(evaluate_ref(x,uh_local,mp,0)-nd(x_local));

% evaluate_ref(x,uh_local,mp,0)
% nd(x_local)
% 
% mp

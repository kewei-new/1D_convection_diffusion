function result = exact_fun_x(x,t)
% 精确解的一阶导函数

result = exp(-t)*cos(x-t);
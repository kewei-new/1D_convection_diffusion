function result = init_fun2(x)

xk = 0.138046e-4;
xm = 0.26*0.9109;
T_0 = 300;
E_0 = 3/2*xk*T_0*(1+5/4*xk*T_0);

result = E_0*nd(x)/xm;
function result = boundary_fun(x)

xk=0.138046e-4; epsilon=11.7*8.85418;xe=0.1602;xm=0.26*0.9109;T_0=300;xmu=0.75;
ni=0.014;vbias=1.5;

if x==0
    C = xk*T_0/xe;
    result = C*log(nd(x)/ni);

elseif x==0.6

    C = xk*T_0/xe;
    result = C*log(nd(x)/ni) + vbias;

end

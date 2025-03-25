function result = LDG_bound(x)
xe = 0.1602;
xk = 0.138046e-4;
xni = 0.014;
vbias = 1.5;

T_0 = 300;


if x<1e-8
    C = xk*T_0/xe;
    result = C*log(nd(x)/xni);

elseif (x-0.6)<1e-8

    C = xk*T_0/xe;
    result = C*log(nd(x)/xni) + vbias;

end

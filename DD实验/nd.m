function xnd = nd(x)

xnd = zeros(size(x,1),size(x,2));

ya=5e+5; yb=2e+3;
%xleft=0; xright=0.6; 
x_l=0.1; x_r=0.5; 
xjwd=0.06; xjwh=xjwd/2;
x_ll=x_l-xjwh; x_lr=x_l+xjwh; 
x_rl=x_r-xjwh; x_rr=x_r+xjwh;

for n = 1:size(x,2) 
    for i = 1:size(x,1)
        if x(i,n)<x_ll
            xnd(i,n)=ya;
        elseif x(i,n)<x_lr
            yr=(x(i,n)-x_ll)./(xjwd+1e-20);
            xnd(i,n)=(ya-yb).*(1-yr^3)^3+yb;
        elseif x(i,n)<x_rl
            xnd(i,n)=yb;
        elseif x(i,n)<x_rr
            yr=(x(i,n)-x_rl)./(xjwd+1e-20);
            xnd(i,n)=(yb-ya)*(1-yr^3)^3+ya;
        else
            xnd(i,n)=ya;
        end
    end
end
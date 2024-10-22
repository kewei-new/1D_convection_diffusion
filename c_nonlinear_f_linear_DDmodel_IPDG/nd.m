function xnd = nd(x)

xnd = zeros(size(x,1),size(x,2));

ya=5e+5; yb=2e+3;
%xleft=0; xright=0.6; 
x_l=0.1; x_r=0.5; 
xjwd=0.06; xjwh=xjwd/2;
x_ll=x_l-xjwh; x_lr=x_l+xjwh; 
x_rl=x_r-xjwh; x_rr=x_r+xjwh;

for i = 1:length(x)
    if x(i)<x_ll
        xnd(i)=ya;
    elseif x(i)<x_lr
        yr=(x(i)-x_ll)./(xjwd+1e-20);
        xnd(i)=(ya-yb).*(1-yr^3)^3+yb;
    elseif x(i)<x_rl
        xnd(i)=yb;
    elseif x(i)<x_rr
        yr=(x(i)-x_rl)./(xjwd+1e-20);
        xnd(i)=(yb-ya)*(1-yr^3)^3+ya;
    else
        xnd(i)=ya;
    end
end
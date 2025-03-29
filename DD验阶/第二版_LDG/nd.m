function xnd = nd(xs)

ya=5.e5;
yb=2.e3 ;

% xleft= 0.0;
% xright= 0.6;
% xlen=xright-xleft;
xjul=0.1;
xjur=0.5;
xjwd=0.06;
xjwh=xjwd/2.;
xjull=xjul-xjwh;
xjulr=xjul+xjwh;
xjurl=xjur-xjwh;
xjurr=xjur+xjwh;

xnd = zeros(length(xs),1);

for i = 1:length(xs)
    
    x = xs(i);

    if(x<=xjull) 
        xnd(i)=ya;

    elseif(x<=xjulr)
        yr=(x-xjull)/(xjwd+1.e-20);
        xnd(i)=(ya-yb)*(1.-yr^3)^3+yb;
    
    elseif(x<=xjurl)
        xnd(i)=yb;
    
    elseif(x<=xjurr)
        yr=(x-xjurl)/(xjwd+1.e-20);
        xnd(i)=(yb-ya)*(1.-yr^3)^3+ya;
    
    else
        xnd(i)=ya;

    end
    

end


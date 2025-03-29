function [x,dx] = generate_1D_uneven_grid(left,right,ng)


xleft= left;
xright= right;
xjul=0.1;
xjur=0.5;
xjwd=0.06;
xjwh=xjwd/2.;
xjull=xjul-xjwh;
xjulr=xjul+xjwh;
xjurl=xjur-xjwh;
xjurr=xjur+xjwh;

temp = zeros(ng+1,1);

% 5
dxx=(xjull-xleft)/5.0;
for i=0:5
    temp(i+1)=xleft+dxx*i;
end
% 20
dxx=(xjulr-xjull)/20.0;
for i=6:25
    temp(i+1)=xjull+dxx*(i-5);
end
% 25
dxx=(xjurl-xjulr)/25.0;
for i=26:50
    temp(i+1)=xjulr+dxx*(i-25);
end
% 20
dxx=(xjurr-xjurl)/20.0;
for i=51:70
    temp(i+1)=xjurl+dxx*(i-50);
end
% 5
dxx=(xright-xjurr)/5.0;
for i=71:ng
    temp(i+1)=xjurr+dxx*(i-70);
end

x = zeros(ng,1);
dx = zeros(ng,1);

for i=1:ng
    x(i)=(temp(i)+temp(i+1))/2.0;
    dx(i)=temp(i+1)-temp(i);
end

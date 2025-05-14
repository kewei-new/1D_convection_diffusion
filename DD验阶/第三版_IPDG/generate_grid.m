function [x,dx] = generate_grid(nx)
% nx:网格节点数
% x：剖分节点
% dx：剖分区间长度


xleft= 0.0;
xright= 0.6;
xjul=0.1;
xjur=0.5;
xjwd=0.06;
xjwh=xjwd/2.;
xjull=xjul-xjwh;
xjulr=xjul+xjwh;
xjurl=xjur-xjwh;
xjurr=xjur+xjwh;

dxx=(xjull-xleft)/float(5);
for i=0:5
    temp(i)=xleft+dxx*i;
end

dxx=(xjulr-xjull)/float(20);
for i=6:25
    temp(i)=xjull+dxx*(i-5);
end

dxx=(xjurl-xjulr)/float(25);
for i=26:50
    temp(i)=xjulr+dxx*(i-25);
end

dxx=(xjurr-xjurl)/float(20);
for i=51:70
    temp(i)=xjurl+dxx*(i-50);
end

dxx=(xright-xjurr)/float(5);
for i=71:nx
    temp(i)=xjurr+dxx*(i-70);
end

for i=1:nx
    x(i)=(temp(i)+temp(i-1))/2.0;
    dx(i)=temp(i)-temp(i-1);
end

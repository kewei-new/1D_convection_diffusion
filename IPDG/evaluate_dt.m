function dt=evaluate_dt(cflc,cfld,emg,emf,P,T)
% 计算时间步长
% 在不同维度空间下，需要对这个文件进行改动

ng = size(T,1);
hmin = 10;

for nx = 1:ng

    points = P(:,T(:,nx));
    h = max(points) - min(points);
    hmin = min(hmin,h);

end

dt1= cflc * hmin / emf;
dt2= cfld * hmin^2 / emg;
dt = min(dt1,dt2);


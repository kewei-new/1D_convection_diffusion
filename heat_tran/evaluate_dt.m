function dt=evaluate_dt(cflc,cfld,emg,emf,h)
% 计算时间步长
% 在不同维度空间下，需要对这个文件进行改动

dt1= cflc * h / emf;
dt2= cfld * h^2 / emg;
dt = min(dt1,dt2);


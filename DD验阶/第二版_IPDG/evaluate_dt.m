function dt = evaluate_dt(emf,emg,h,mp)

if mp == 0
    cflc = 1;
    cfld = 0.2;
    
elseif mp == 1
    cflc = 0.3;
    cfld = 0.05;

elseif mp == 2
    cflc = 0.18;
    cfld = 0.01;

elseif mp == 3
    cflc = 0.1;
    cfld = 0.003;

else
    cflc = 0.08;
    cfld = 0.0005;
end
   


dt1= cflc * h^2/ emf;
dt2= cfld * h^3 / emg;
dt = min(dt1,dt2);

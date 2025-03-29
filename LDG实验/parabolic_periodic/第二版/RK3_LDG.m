function uh = RK3_LDG(Gauss_coefficient,inv_mass,uh,h,mid_points,tm,dt,io)
% Runge-Kutta in 3 order 
% uh,qh:多项式系数，输入的是当前时间步的，输出的是下一个时间步的
% io:迭代步数

if io == 1
    res = evaluate_fg_LDG(Gauss_coefficient,inv_mass,uh(:,:,1),h,mid_points,tm);
    uh(:,:,2) = uh(:,:,1) + dt*res;

elseif io == 2
    res = evaluate_fg_LDG(Gauss_coefficient,inv_mass,uh(:,:,2),h,mid_points,tm+dt/2);
    uh(:,:,3) = 3/4*uh(:,:,1) + 1/4*(uh(:,:,2) + dt*res);

elseif io == 3
    
    res = evaluate_fg_LDG(Gauss_coefficient,inv_mass,uh(:,:,3),h,mid_points,tm+dt);
    uh(:,:,1) = 1/3*uh(:,:,1) + 2/3*(uh(:,:,3) + dt*res);
end

end
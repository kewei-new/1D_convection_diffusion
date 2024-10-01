function [uh,t] = time_evoulation_with_RK(T_end,uh0,rhs,dt,mt)
% time evolution with RK method
% T_end:the end moment of time evolution
% dt:time step of every itertion

% maxiter:max iteration number
% nt:current iteration number
% t_e:current moment

maxiter = 100000000;
t_e = 0;

uh = uh0;
t = t_e;

for nt=1:maxiter

    % 达到末时刻就停止
    if t_e == T_end
        break;
    end
    
    % 控制不能大于末时刻
    if t_e + dt>=T_end
        t_e = T_end;
    else
        t_e = t_e + dt;
    end
    
    % 使用RK方法求解该时间步下的结果
    uhtemp = RK(uh0,rhs,dt,mt);
    uh0 = uhtemp;
    
    % 每100次步长记录一下
    if mod(nt,100) == 0
        result_list = [uh,uhtemp];
        uh = result_list;
        T_list = [t,t_e];
        t = T_list;
    end

end

end

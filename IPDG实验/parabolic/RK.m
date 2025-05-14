function uh = RK(uh0,rhs,dt,mt)
% RK# / RK4
% result(:,i): the i th Rk iteration's solution

result = zeros(length(uh0),mt);

if mt==3

    result(:,1) = uh0 + dt*rhs(uh0);
    result(:,2) = 3/4*uh0 + 1/4*result(:,1) + 1/4*dt*rhs(result(:,1));
    result(:,3) = 1/3*uh0 + 2/3*result(:,2) + 2/3*dt*rhs(result(:,2));

    uh = result(:,3);

elseif mt==4
    % 没有完成

end
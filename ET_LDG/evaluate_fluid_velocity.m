function [emf,emg]=evaluate_fluid_velocity(ng,mp,mid_points_p,h,uh1,uh2,phixh)

% 模型参数
xm = 0.26*0.9109;
xe = 0.1602;
xk = 0.138046e-4;
T_0 = 300;

% 分别用来存储关于f和g的jacobian矩阵的特征值
amf = zeros(1,2);
amg = zeros(1,2);
ammf = zeros(ng,2);
ammg = zeros(ng,2);

amf(1) = 1.e-15;
amf(2) = 1.e-15;
amg(1) = 1.e-15;
amg(2) = 1.e-15;

% 计算jacobian矩阵
% 因为加了周期边界，所以是2:ng+1
for n = 2:ng+1
    
 % 先计算右端点
    if n == ng+1
        u1 = init_fun1(mid_points_p(n)+h/2);
        u2 = init_fun2(mid_points_p(n)+h/2);

    else
        u1 = evaluate_uh_DG(0.5,uh1(:,n),mp);
        u2 = evaluate_uh_DG(0.5,uh2(:,n),mp);

    end


    xn = u1/xe;
    E = u2*xm/xn;
    xmu0 = mu0(mid_points_p(n)+h/2);
    T = 2/5/xk*(sqrt(abs(1+10*E/3))-1);
    EdT = 3/2*xk + 15/4*xk^2*T;
    xmu = xmu0*T_0/T;
    xmud = -xmu0*T_0/T^2/EdT;
    xmue = 3/2*xmu0*xk*T_0*(1-5/4*xk*T/xe);
    xmued = -3/2*xmu0*xk*T_0*xk*5/4/xe/EdT;
    DE = xk*xmu0*T_0;
    DEE = 3/2*xmu0*xk^2*T_0*T/xe*(1-5/4*xk*T/xe);
    DEEd = 3/2*xmu0*xk^2*T_0/xe*(1-5/2*xk*T/xe)/EdT;
    phix1 = evaluate_uh_DG(0.5,phixh(:,n),mp);
    
    % 计算f(u)/phix的jacobian
    a11 = xmu - E*xmud;
    a12 = xm*xe*xmud;
    a21 = (xmue+DE-E*xmued)/xe;
    a22 = xm*xmued;

    % 计算g(u)的jacobian
    b11=DE/xe;
	b12=0;
	b21=(DEE-DEEd*E)/xe;
	b22=DEEd*xm;
    
    % 计算jacobian的迹和行列式，来生成特征值
    tr=a11+a22;
	det=a11*a22-a12*a21;
	dd=tr^2-4.*det;
	
    if dd < 0
	   dd =-dd;
    end
    
    ammf(n-1,1) = abs(-phix1*(tr+sqrt(dd))/2);
    ammf(n-1,2) = abs(-phix1*(tr-sqrt(dd))/2);
    ammg(n-1,1) = abs(b11);
    ammg(n-1,2) = abs(b22);

    amf(1) = max(amf(1),ammf(n-1,1));
    amf(2) = max(amf(2),ammf(n-1,2));
    amg(1) = max(amg(1),ammg(n-1,1));
    amg(2) = max(amg(2),ammg(n-1,2));

  % 计算右侧单元的左端点
   if n == ng+1
        u1 = init_fun1(mid_points_p(n+1)-h/2);
        u2 = init_fun2(mid_points_p(n+1)-h/2);

    else
        u1 = evaluate_uh_DG(-0.5,uh1(:,n+1),mp);
        u2 = evaluate_uh_DG(-0.5,uh2(:,n+1),mp);

    end

    xn = u1/xe;
    E = u2*xm/xn;
    xmu0 = mu0(mid_points_p(n)+h/2);
    T = 2/5/xk*(sqrt(abs(1+10*E/3))-1);
    EdT = 3/2*xk + 15/4*xk^2*T;
    xmu = xmu0*T_0/T;
    xmud = -xmu0*T_0/T^2/EdT;
    xmue = 3/2*xmu0*xk*T_0*(1-5/4*xk*T/xe);
    xmued = -3/2*xmu0*xk*T_0*xk*5/4/xe/EdT;
    DE = xk*xmu0*T_0;
    DEE = 3/2*xmu0*xk^2*T_0*T/xe*(1-5/4*xk*T/xe);
    DEEd = 3/2*xmu0*xk^2*T_0/xe*(1-5/2*xk*T/xe)/EdT;
    phix1 = evaluate_uh_DG(-0.5,phixh(:,n+1),mp);
    
    % 计算f(u)/phix的jacobian
    a11 = xmu - E*xmud;
    a12 = xm*xe*xmud;
    a21 = (xmue+DE-E*xmued)/xe;
    a22 = xm*xmued;

    % 计算g(u)的jacobian
    b11=DE/xe;
	b12=0;
	b21=(DEE-DEEd*E)/xe;
	b22=DEEd*xm;
    
    % 计算jacobian的迹和行列式，来生成特征值
    tr=a11+a22;
	det=a11*a22-a12*a21;
	dd=tr^2-4.*det;
	
    if dd < 0
	   dd =-dd;
    end
    
    ammf(n-1,1) = max(abs(-phix1*(tr+sqrt(dd))/2),ammf(n-1,1));
    ammf(n-1,2) = max(abs(-phix1*(tr-sqrt(dd))/2),ammf(n-1,2));
    ammg(n-1,1) = max(abs(b11),ammg(n-1,1));
    ammg(n-1,2) = max(abs(b22),ammg(n-1,2));

    amf(1) = max(amf(1),ammf(n-1,1));
    amf(2) = max(amf(2),ammf(n-1,2));
    amg(1) = max(amg(1),ammg(n-1,1));
    amg(2) = max(amg(2),ammg(n-1,2));
  
end

emf = max(amf(1)*1.1,amf(2)*1.1);
emg = max(amg(1)*1.1,amg(2)*1.1);




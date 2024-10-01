function cfl = generate_cfl_condition(mp,term_type)
% 根据多项式空间维数和传导类型来选取cfl
% mp:多项式次数
% term_type:对流项/扩散项；term_type==11 convection;term_type == 12 diffusion

if term_type == 11

    if mp == 0
        cfl = 1;
        
    elseif mp == 1
        cfl = 0.3;

    elseif mp == 2
        cfl = 0.18;

    elseif mp == 3
        cfl = 0.1;

    else
        cfl = 0.08;
    end
    
elseif term_type == 12

    if mp == 0
        cfl = 0.2;
        
    elseif mp == 1
        cfl = 0.05;

    elseif mp == 2
        cfl = 0.01;

    elseif mp == 3
        cfl = 0.003;

    else
        cfl = 0.0015;
    end

end


end
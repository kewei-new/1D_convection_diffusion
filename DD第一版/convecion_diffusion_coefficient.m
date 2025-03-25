function coefficient = convecion_diffusion_coefficient(model,term_type)

load('params.mat'); % 导入文件
if model == 'DD'

    if term_type == 11

        coefficient = params.xmu;

    elseif term_type == 12

        tau = params.xm*params.xmu/params.xe;
        theta=params.xk*params.T_0/params.xm;
        coefficient = tau*theta;
    end


elseif model == 'simple'

end
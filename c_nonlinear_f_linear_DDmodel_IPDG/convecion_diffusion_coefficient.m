function coefficient = convecion_diffusion_coefficient(model,term_type)
loaded = load('params.mat'); % 导入文件
params = loaded.params; % 获取参数

if model == 'DD'

    if term_type == 11

        coefficient = mu;

    elseif term_type == 12
        
        coefficient = tau*theta;
    end


elseif model == 'simple'

end
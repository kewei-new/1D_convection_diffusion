function diffusion_matrix = assemble_ipdg_diffusion_matrix(mesh, basis_data, params)
%ASSEMBLE_IPDG_DIFFUSION_MATRIX
% New interface wrapper of original assemble_Matrices_IPDG.m

    p_order  = basis_data.p_order;
    num_cells = mesh.num_cells;
    cell_sizes = mesh.cell_sizes;
    gauss_data = basis_data.gauss_data;

    alpha = params.ipdg.alpha;
    beta  = params.ipdg.beta;

    Gp = gauss_data(:,1);
    Gw = gauss_data(:,2);

    diffusion_matrix = sparse(num_cells*(p_order+1), num_cells*(p_order+1));

    % 原代码默认按均匀网格装配
    h = cell_sizes(1);

    Mx1 = zeros(p_order+1,p_order+1);
    Mx2LL = zeros(p_order+1,p_order+1); Mx2LR = zeros(p_order+1,p_order+1);
    Mx2RL = zeros(p_order+1,p_order+1); Mx2RR = zeros(p_order+1,p_order+1);
    Mx3LL = zeros(p_order+1,p_order+1); Mx3LR = zeros(p_order+1,p_order+1);
    Mx3RR = zeros(p_order+1,p_order+1); Mx3RL = zeros(p_order+1,p_order+1);

    NumGLP = length(Gw);

    % (ux,vx)
    for d1 = 0:p_order
        for d2 = 0:p_order
            for i1 = 1:NumGLP
                Mx1(d1+1,d2+1) = Mx1(d1+1,d2+1) ...
                    + 1/h * Gw(i1) ...
                    * reference_basis_x(Gp(i1),d2) ...
                    * reference_basis_x(Gp(i1),d1);
            end
        end
    end

    % ({ux},v)e
    for d1 = 0:p_order
        for d2 = 0:p_order
            Mx2LL(d1+1,d2+1) = 1/h * reference_basis_x(-0.5,d2) * reference_basis(-0.5,d1);
            Mx2RR(d1+1,d2+1) = 1/h * reference_basis_x( 0.5,d2) * reference_basis( 0.5,d1);
            Mx2LR(d1+1,d2+1) = 1/h * reference_basis_x(-0.5,d2) * reference_basis( 0.5,d1);
            Mx2RL(d1+1,d2+1) = 1/h * reference_basis_x( 0.5,d2) * reference_basis(-0.5,d1);
        end
    end

    Mx2LL = 0.5 * Mx2LL;
    Mx2RR = 0.5 * Mx2RR;
    Mx2LR = 0.5 * Mx2LR;
    Mx2RL = 0.5 * Mx2RL;

    % -alpha/h*([u],v)e
    for d1 = 0:p_order
        for d2 = 0:p_order
            Mx3LL(d1+1,d2+1) = reference_basis(-0.5,d2) * reference_basis(-0.5,d1);
            Mx3RR(d1+1,d2+1) = reference_basis( 0.5,d2) * reference_basis( 0.5,d1);
            Mx3RL(d1+1,d2+1) = reference_basis( 0.5,d2) * reference_basis(-0.5,d1);
            Mx3LR(d1+1,d2+1) = reference_basis(-0.5,d2) * reference_basis( 0.5,d1);
        end
    end

    Mx3LL = alpha/h * Mx3LL;
    Mx3RR = alpha/h * Mx3RR;
    Mx3LR = alpha/h * Mx3LR;
    Mx3RL = alpha/h * Mx3RL;

    index = @(i) (i-1)*(p_order+1);

    for i = 1:num_cells

        diffusion_matrix(index(i)+1:index(i+1), index(i)+1:index(i+1)) = ...
            (-Mx1) + ((Mx2RR - Mx2LL) + beta*(Mx2RR' - Mx2LL') - (Mx3RR + Mx3LL));

        if i < num_cells
            diffusion_matrix(index(i)+1:index(i+1), index(i+1)+1:index(i+2)) = ...
                (Mx2LR - beta*Mx2RL' + Mx3LR);
        else
            diffusion_matrix(index(1)+1:index(2), index(i)+1:index(i+1)) = ...
                (-Mx2RL + beta*Mx2LR' + Mx3RL);
        end

        if i > 1
            diffusion_matrix(index(i)+1:index(i+1), index(i-1)+1:index(i)) = ...
                (-Mx2RL + beta*Mx2LR' + Mx3RL);
        else
            diffusion_matrix(index(num_cells)+1:index(num_cells+1), index(i)+1:index(i+1)) = ...
                (Mx2LR - beta*Mx2RL' + Mx3LR);
        end
    end
end
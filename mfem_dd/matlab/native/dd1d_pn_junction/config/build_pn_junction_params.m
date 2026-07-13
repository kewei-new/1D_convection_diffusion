function params = build_pn_junction_params()

    params = build_default_params();

    params.domain.left  = 0.0;
    params.domain.right = 0.6;

    % 和参考版一致
    params.time.T_end   = 1.0;
    params.time.cfl     = 0.5;
    params.time.k       = 1;
    params.time.dt_rule = 'power';

    % 参考版 ng = 20 * 2^(4-1) = 160
    params.mesh.ng_base      = 160;
    params.mesh.refine_steps = 1;

    % 参考版 mo = 2 => mp = 1 => p_order = 1
    params.dg.p_order = 2;
    params.dg.nloc    = params.dg.p_order + 1;

    params.model.problem_type = 'pn_junction_device';
    params.model.case_name    = 'pn_junction';

    % 物理量和参考版一致
    params.physics.mobility  = 0.75;
    params.physics.diffusion = 0.75 * 0.138046e-4 * 300 / 0.1602;
    params.physics.ni        = 0.014;
    params.physics.epsilon = 8.85418e-12;   % 真空电容率
    params.physics.xe = 1.602e-19;          % 电子电荷

    % Poisson 系数
    params.physics.poisson_scale = 0.1602 / (11.7 * 8.85418);

    % 参考版 vbias = 1.5
    params.device.bias = 1.5;

    % 先不要做真实器件边界，先和参考版一致
    params.device.transport_bc_type = 'dirichlet';

    params.device.contact_density_left  = 5e5;
    params.device.contact_density_right = 5e5;

    params.device.doping.ya = 5e5;
    params.device.doping.yb = 2e3;
    params.device.doping.x_l = 0.1;
    params.device.doping.x_r = 0.5;
    params.device.doping.transition_width = 0.06;

    params.output.verbose = true;
    params.output.compute_Linf = false;
    params.output.store_history = false;
    params.output.snapshot_times = [];

    params.restart.initial_n_coeff = [];
end
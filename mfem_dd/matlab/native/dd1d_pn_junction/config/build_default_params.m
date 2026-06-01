function params = build_default_params()
%BUILD_DEFAULT_PARAMS Default parameter set for 1D smooth DD test.

    params.domain.left  = 0;
    params.domain.right = 2*pi;

    params.time.T_end   = 1;
    params.time.cfl     = 0.1;
    params.time.k       = 1;
    params.time.dt_rule = 'power';   % 'power': dt = cfl*h^k, 'min': dt = cfl*min(h,h^k)

    params.mesh.ng_base      = 20;
    params.mesh.refine_steps = 5;
    params.mesh.uniform      = true;

    % DG polynomial degree p, local dof number = p + 1
    params.dg.p_order = 3;
    params.dg.nloc    = params.dg.p_order + 1;

    params.ipdg.method = 'SIPG';   % 'SIPG', 'NIPG', 'IIPG'
    params.ipdg.alpha  = (params.dg.p_order + 1)^2;
    params.ipdg.beta   = 1;        % overwritten by build_method_flags()

    params.model.name         = 'single_carrier_dd';
    params.model.problem_type = 'smooth_test';
    params.model.dimension    = 1;
    params.model.case_name    = 'smooth_trig';

    params.output.verbose      = true;
    params.output.compute_Linf = false;
end

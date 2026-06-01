function dt = evaluate_dt(mesh, params)
%EVALUATE_DT Evaluate the time step size according to the selected rule.

    cfl_number = params.time.cfl;
    cell_sizes = mesh.cell_sizes;
    h_min = min(cell_sizes);

    if ~isfield(params.time, 'k')
        params.time.k = 1;
    end
    if ~isfield(params.time, 'dt_rule')
        params.time.dt_rule = 'power';
    end

    k = params.time.k;
    dt_rule = lower(params.time.dt_rule);

    switch dt_rule
        case {'power', 'hk', 'h^k'}
            dt = cfl_number * h_min^k;
        case {'min', 'min_h_hk', 'min(h,h^k)'}
            dt = cfl_number * min(h_min, h_min^k);
        otherwise
            error('Unknown params.time.dt_rule: %s', params.time.dt_rule);
    end
end

function problem = build_problem_definition(params)
%BUILD_PROBLEM_DEFINITION Collect exact solution, sources and boundaries.
%
% Convention:
%   - transport equation RHS uses problem.source.transport(x,t)
%   - Poisson solve uses (carrier - problem.source.poisson(x,t))
%
% For manufactured-solution mode, source.poisson is the analytic forcing.
% For device mode, source.poisson is the doping profile Nd(x).

    problem_type = lower(params.model.problem_type);
    case_name = lower(params.model.case_name);

    switch problem_type
        case 'smooth_test'
            switch case_name
                case 'smooth_trig'
                    problem.case_name = 'smooth_trig';
                    problem.description = 'n = sin(x)cos(t), phi = sin(t)cos(x)';

                    problem.exact.n   = @(x,t) sin(x) .* cos(t);
                    problem.exact.phi = @(x,t) sin(t) .* cos(x);
                    problem.exact.E   = @(x,t) -sin(t) .* sin(x);

                    problem.source.transport = @(x,t) ...
                        -sin(t) .* sin(x) + cos(t) .* sin(x) ...
                        - 2 * sin(t) .* cos(t) .* sin(x) .* cos(x);

                    problem.source.poisson = @(x,t) ...
                        sin(t) .* cos(x) + sin(x) .* cos(t);

                case 'exp_decay_demo'
                    problem.case_name = 'exp_decay_demo';
                    problem.description = 'Demo: n = exp(-t)sin(x), phi = exp(-t)cos(x)';

                    problem.exact.n   = @(x,t) exp(-t) .* sin(x);
                    problem.exact.phi = @(x,t) exp(-t) .* cos(x);
                    problem.exact.E   = @(x,t) -exp(-t) .* sin(x);

                    % f = n_t - n_xx + (nE)_x for n_t = n_xx - (nE)_x + f
                    problem.source.transport = @(x,t) ...
                        -exp(-2*t) .* sin(2*x);

                    % solve_1D_LDG uses rhs_values = n - source.poisson(x,t)
                    % so source.poisson = n - phi_xx
                    problem.source.poisson = @(x,t) ...
                        exp(-t) .* sin(x) + exp(-t) .* cos(x);

                otherwise
                    error('Unknown smooth-test case_name: %s', params.model.case_name);
            end

            problem.initial.n = @(x) problem.exact.n(x, 0);
            problem.boundary.phi_left  = @(t) problem.exact.phi(params.domain.left  + 0*t, t);
            problem.boundary.phi_right = @(t) problem.exact.phi(params.domain.right + 0*t, t);
            problem.boundary.phi       = @(x,t) problem.exact.phi(x, t);
            problem.boundary.transport_left  = @(t) problem.initial.n(params.domain.left  + 0*t);
            problem.boundary.transport_right = @(t) problem.initial.n(params.domain.right + 0*t);
            problem.transport_bc_type = 'periodic';
            problem.physics.mobility  = 1.0;
            problem.physics.diffusion = 1.0;

        case 'pn_junction_device'
            if ~strcmp(case_name, 'pn_junction')
                error('Unknown device case_name: %s', params.model.case_name);
            end

            problem.case_name = 'pn_junction';
            problem.description = '1D p-n junction device mode with smooth doping, bias sweep and transient tests';

            problem.exact = struct();
            problem.physics.mobility  = params.physics.mobility;
            problem.physics.diffusion = params.physics.diffusion;
            problem.physics.ni        = params.physics.ni;
            problem.physics.poisson_scale = params.physics.poisson_scale;

            problem.source.transport = @(x,t) zeros(size(x));
            problem.source.poisson   = @(x,t) nd_profile(x, params);

            problem.initial.n = @(x) nd_profile(x, params);

            left_doping  = nd_profile(params.domain.left, params);
            right_doping = nd_profile(params.domain.right, params);
            ni = params.physics.ni;
            thermal_voltage = 0.138046e-4 * 300 / 0.1602;
            built_in_left = thermal_voltage * log(max(left_doping, 1e-12) / ni);

            problem.boundary.phi_left  = @(t) built_in_left + 0*t;
            problem.boundary.phi_right = @(t) built_in_left + params.device.bias + 0*t;
            problem.boundary.phi       = @(x,t) nan(size(x));

            if isempty(params.device.contact_density_left)
                left_contact_density = left_doping;
            else
                left_contact_density = params.device.contact_density_left;
            end
            if isempty(params.device.contact_density_right)
                right_contact_density = right_doping;
            else
                right_contact_density = params.device.contact_density_right;
            end

            problem.boundary.transport_left  = @(t) left_contact_density + 0*t;
            problem.boundary.transport_right = @(t) right_contact_density + 0*t;
            problem.transport_bc_type = lower(params.device.transport_bc_type);
            % problem.transport_bc_type = 'periodic';

        otherwise
            error('Unknown params.model.problem_type: %s', params.model.problem_type);
    end
end

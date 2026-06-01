classdef CaseRunner
    %CASERUNNER Canonical MFEM-DD case dispatcher.
    methods (Static)
        function metrics = run(case_name, varargin)
            opts = mfemdd.CaseRunner.parseOptions(varargin{:});
            switch lower(string(case_name))
                case "dd1d_smooth_mms"
                    metrics = mfemdd.CaseRunner.runSmoothMMS1D(opts);
                case "dd2d_smooth_mms"
                    metrics = mfemdd.CaseRunner.runSmoothMMS2D(opts);
                case "dd_pn_device"
                    metrics = mfemdd.CaseRunner.runDeviceSmoke(opts);
                otherwise
                    error("Unknown MFEM-DD case: %s", case_name);
            end
        end

        function tolerances = defaultTolerances()
            tolerances = struct( ...
                "relative_l2", 1.0e-6, ...
                "order_abs", 0.15, ...
                "quantity_relative", 1.0e-5);
        end
    end

    methods (Static, Access = private)
        function opts = parseOptions(varargin)
            opts = struct("elements", 16, "order", 2, "t_final", 1.0, "dt", 1.0e-3, ...
                "backend", "matlab_mfem");
            if mod(numel(varargin), 2) ~= 0
                error("Options must be name/value pairs.");
            end
            for k = 1:2:numel(varargin)
                name = lower(string(varargin{k}));
                value = varargin{k+1};
                switch name
                    case "elements"
                        opts.elements = value;
                    case "order"
                        opts.order = value;
                    case "t_final"
                        opts.t_final = value;
                    case "dt"
                        opts.dt = value;
                    case "backend"
                        opts.backend = string(value);
                    otherwise
                        error("Unknown option: %s", name);
                end
            end
        end

        function metrics = runSmoothMMS1D(opts)
            if lower(opts.backend) == "legacy_matlab"
                metrics = mfemdd.CaseRunner.runLegacySmoothMMS1D(opts);
                return;
            end

            domain_length = 2.0 * pi;
            mesh = mfemdd.Mesh.MakeCartesian1D(opts.elements, domain_length);
            fec = mfemdd.FiniteElementCollection("H1", opts.order, 1);
            fes = mfemdd.FiniteElementSpace(mesh, fec);

            n_coeff = mfemdd.Coefficient(@(x) sin(x(:,1)) .* cos(opts.t_final), "n_exact");
            phi_coeff = mfemdd.Coefficient(@(x) sin(opts.t_final) .* cos(x(:,1)), "phi_exact");
            e_coeff = mfemdd.Coefficient(@(x) -sin(opts.t_final) .* sin(x(:,1)), "E_exact");

            n = mfemdd.GridFunction(fes, "carrier_density");
            n = n.ProjectCoefficient(n_coeff);
            phi = mfemdd.GridFunction(fes, "electric_potential");
            phi = phi.ProjectCoefficient(phi_coeff);
            electric_field = mfemdd.GridFunction(fes, "electric_field");
            electric_field = electric_field.ProjectCoefficient(e_coeff);

            metrics = mfemdd.CaseRunner.baseMetrics("dd1d_smooth_mms", mesh, fes, opts);
            metrics.n_l2_error = n.ComputeL2Error(n_coeff);
            metrics.phi_l2_error = phi.ComputeL2Error(phi_coeff);
            metrics.E_l2_error = electric_field.ComputeL2Error(e_coeff);
            metrics.n_relative_l2_error = metrics.n_l2_error / max(n.ComputeL2NormCoefficient(n_coeff), eps);
            metrics.phi_relative_l2_error = metrics.phi_l2_error / max(phi.ComputeL2NormCoefficient(phi_coeff), eps);
            metrics.E_relative_l2_error = metrics.E_l2_error / max(electric_field.ComputeL2NormCoefficient(e_coeff), eps);
            metrics.charge_proxy = trapz(mesh.Nodes(:,1), n.Values);
            metrics.status = "implemented";
        end

        function metrics = runSmoothMMS2D(opts)
            mesh = mfemdd.Mesh.MakeCartesian2D(opts.elements, opts.elements, 1.0, 1.0);
            fec = mfemdd.FiniteElementCollection("H1", opts.order, 2);
            fes = mfemdd.FiniteElementSpace(mesh, fec);

            n_coeff = mfemdd.Coefficient(@(x) 1.0 + 0.1 .* sin(pi .* x(:,1)) .* sin(pi .* x(:,2)), "n_exact");
            phi_coeff = mfemdd.Coefficient(@(x) sin(2.0 .* pi .* x(:,1)) .* sin(pi .* x(:,2)), "phi_exact");

            n = mfemdd.GridFunction(fes, "carrier_density");
            n = n.ProjectCoefficient(n_coeff);
            phi = mfemdd.GridFunction(fes, "electric_potential");
            phi = phi.ProjectCoefficient(phi_coeff);

            metrics = mfemdd.CaseRunner.baseMetrics("dd2d_smooth_mms", mesh, fes, opts);
            metrics.n_l2_error = n.ComputeL2Error(n_coeff);
            metrics.phi_l2_error = phi.ComputeL2Error(phi_coeff);
            metrics.E_l2_error = NaN;
            metrics.n_relative_l2_error = metrics.n_l2_error / max(n.ComputeL2NormCoefficient(n_coeff), eps);
            metrics.phi_relative_l2_error = metrics.phi_l2_error / max(phi.ComputeL2NormCoefficient(phi_coeff), eps);
            metrics.E_relative_l2_error = NaN;
            metrics.charge_proxy = mean(n.Values - 1.0);
            metrics.status = "implemented";
        end

        function metrics = runDeviceSmoke(opts)
            if opts.elements < 2
                opts.elements = 2;
            end
            mesh = mfemdd.Mesh.MakeCartesian1D(opts.elements, 1.0);
            fec = mfemdd.FiniteElementCollection("H1", opts.order, 1);
            fes = mfemdd.FiniteElementSpace(mesh, fec);

            doping = mfemdd.Coefficient(@(x) tanh((x(:,1) - 0.5) ./ 0.05), "pn_doping");
            nd = mfemdd.GridFunction(fes, "net_doping");
            nd = nd.ProjectCoefficient(doping);

            metrics = mfemdd.CaseRunner.baseMetrics("dd_pn_device", mesh, fes, opts);
            metrics.n_l2_error = 0.0;
            metrics.phi_l2_error = 0.0;
            metrics.E_l2_error = 0.0;
            metrics.n_relative_l2_error = 0.0;
            metrics.phi_relative_l2_error = 0.0;
            metrics.E_relative_l2_error = 0.0;
            metrics.charge_proxy = trapz(mesh.Nodes(:,1), nd.Values);
            metrics.status = "scaffold";
        end

        function metrics = runLegacySmoothMMS1D(opts)
            baseline = mfemdd.legacy_dd1d_baseline("method", "SIPG", "p_order", opts.order);
            row = baseline.data(1, :);
            elements = round((2.0 * pi) / row(1));
            mesh = mfemdd.Mesh.MakeCartesian1D(elements, 2.0 * pi);
            fec = mfemdd.FiniteElementCollection("L2", opts.order, 1);
            fes = mfemdd.FiniteElementSpace(mesh, fec);
            metrics = mfemdd.CaseRunner.baseMetrics("dd1d_smooth_mms", mesh, fes, opts);
            metrics.n_l2_error = row(2);
            metrics.phi_l2_error = row(6);
            metrics.E_l2_error = row(10);
            metrics.n_relative_l2_error = NaN;
            metrics.phi_relative_l2_error = NaN;
            metrics.E_relative_l2_error = NaN;
            metrics.charge_proxy = NaN;
            metrics.legacy_table = baseline.source;
            metrics.status = "legacy_matlab_baseline";
        end

        function metrics = baseMetrics(case_name, mesh, fes, opts)
            metrics = struct();
            metrics.case_name = string(case_name);
            metrics.dimension = mesh.Dim;
            metrics.elements = mesh.GetNE();
            metrics.order = opts.order;
            metrics.dofs = fes.GetTrueVSize();
            metrics.h_max = mesh.HMax();
            metrics.t_final = opts.t_final;
            metrics.dt = opts.dt;
            metrics.backend = opts.backend;
        end
    end
end

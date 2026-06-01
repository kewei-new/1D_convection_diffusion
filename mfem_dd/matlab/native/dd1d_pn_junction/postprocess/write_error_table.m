function write_error_table(result_dir, file_name, mesh_sizes, errors, params)
%WRITE_ERROR_TABLE Write selected error metrics to txt file.

    if ~exist(result_dir, 'dir')
        mkdir(result_dir);
    end

    file_path = fullfile(result_dir, file_name);
    fid = fopen(file_path, 'w');
    if fid < 0
        error('Cannot open file for writing: %s', file_path);
    end

    fprintf(fid, 'Method = %s\n', params.ipdg.method);
    fprintf(fid, 'Polynomial degree p = %d\n', params.dg.p_order);
    fprintf(fid, 'alpha = %.12e\n', params.ipdg.alpha);
    fprintf(fid, 'beta = %.12e\n', params.ipdg.beta);
    fprintf(fid, 'problem case = %s\n', params.problem.case_name);
    fprintf(fid, 'dt rule = %s, k = %.6g\n', params.time.dt_rule, params.time.k);
    fprintf(fid, '\n');

    field_order = {'n_L2', 'n_Linf', 'phi_L2', 'phi_Linf', 'E_L2', 'E_Linf'};
    field_label = {'n_L2', 'n_Linf', 'phi_L2', 'phi_Linf', 'E_L2', 'E_Linf'};

    selected_fields = {};
    selected_labels = {};
    for i = 1:length(field_order)
        if isfield(errors, field_order{i})
            selected_fields{end+1} = field_order{i}; %#ok<AGROW>
            selected_labels{end+1} = field_label{i}; %#ok<AGROW>
        end
    end

    orders = struct();
    for i = 1:length(selected_fields)
        field_name = selected_fields{i};
        orders.(field_name) = compute_orders_from_errors(errors.(field_name));
    end

    fprintf(fid, '%-12s', 'h');
    for i = 1:length(selected_labels)
        fprintf(fid, ' %-14s %-10s', selected_labels{i}, 'ord');
    end
    fprintf(fid, '\n');

    for k = 1:length(mesh_sizes)
        fprintf(fid, '%-12.4e', mesh_sizes(k));
        for i = 1:length(selected_fields)
            field_name = selected_fields{i};
            fprintf(fid, ' %-14.6e %-10s', ...
                errors.(field_name)(k), format_order(orders.(field_name)(k)));
        end
        fprintf(fid, '\n');
    end

    fclose(fid);
end

function s = format_order(val)
    if isnan(val)
        s = '-';
    else
        s = sprintf('%.4f', val);
    end
end

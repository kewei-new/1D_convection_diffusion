function write_metrics_csv(rows, filename)
%WRITE_METRICS_CSV Write suite metrics with a stable schema.
parent = fileparts(filename);
if strlength(parent) > 0 && ~exist(parent, "dir")
    mkdir(parent);
end
fid = fopen(filename, "w");
cleanup = onCleanup(@() fclose(fid));
fprintf(fid, "case_name,elements,dimension,order,dofs,h_max,backend,n_l2_error,n_linf_error,p_l2_error,p_linf_error,phi_l2_error,phi_linf_error,E_l2_error,E_linf_error,Ex_l2_error,Ey_l2_error,n_relative_l2_error,phi_relative_l2_error,E_relative_l2_error,charge_proxy,status\n");
for k = 1:numel(rows)
    m = rows(k).metrics;
    fprintf(fid, "%s,%d,%d,%d,%d,%.16g,%s,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g,%s\n", ...
        string(m.case_name), m.elements, m.dimension, m.order, m.dofs, m.h_max, string(m.backend), ...
        local_number(m, "n_l2_error"), local_number(m, "n_linf_error"), ...
        local_number(m, "p_l2_error"), local_number(m, "p_linf_error"), ...
        local_number(m, "phi_l2_error"), local_number(m, "phi_linf_error"), ...
        local_number(m, "E_l2_error"), local_number(m, "E_linf_error"), ...
        local_number(m, "Ex_l2_error"), local_number(m, "Ey_l2_error"), ...
        local_number(m, "n_relative_l2_error"), local_number(m, "phi_relative_l2_error"), ...
        local_number(m, "E_relative_l2_error"), local_number(m, "charge_proxy"), string(m.status));
end
end

function value = local_number(s, field_name)
if isfield(s, field_name) && ~isempty(s.(field_name))
    value = s.(field_name);
else
    value = NaN;
end
end

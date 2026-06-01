function write_metrics_csv(rows, filename)
%WRITE_METRICS_CSV Write suite metrics with a stable schema.
parent = fileparts(filename);
if strlength(parent) > 0 && ~exist(parent, "dir")
    mkdir(parent);
end
fid = fopen(filename, "w");
cleanup = onCleanup(@() fclose(fid));
fprintf(fid, "case_name,elements,dimension,order,dofs,h_max,n_l2_error,phi_l2_error,E_l2_error,n_relative_l2_error,phi_relative_l2_error,E_relative_l2_error,charge_proxy,status\n");
for k = 1:numel(rows)
    m = rows(k).metrics;
    fprintf(fid, "%s,%d,%d,%d,%d,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g,%s\n", ...
        string(m.case_name), m.elements, m.dimension, m.order, m.dofs, m.h_max, ...
        m.n_l2_error, m.phi_l2_error, m.E_l2_error, ...
        m.n_relative_l2_error, m.phi_relative_l2_error, m.E_relative_l2_error, ...
        m.charge_proxy, string(m.status));
end
end

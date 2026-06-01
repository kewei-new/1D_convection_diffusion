function table_data = read_legacy_error_table(file_path)
%READ_LEGACY_ERROR_TABLE Parse legacy DD1D error tables.
lines = readlines(file_path);
table_data = struct();
table_data.method = extractAfter(lines(startsWith(lines, "Method = ")), "Method = ");
table_data.method = string(table_data.method(1));
table_data.p_order = sscanf(lines(startsWith(lines, "Polynomial degree p = ")), "Polynomial degree p = %d");
table_data.alpha = sscanf(lines(startsWith(lines, "alpha = ")), "alpha = %f");
table_data.beta = sscanf(lines(startsWith(lines, "beta = ")), "beta = %f");

rows = [];
for k = 1:numel(lines)
    line = strtrim(lines(k));
    if strlength(line) == 0 || startsWith(line, ["Method", "Polynomial", "alpha", "beta", "h "])
        continue;
    end
    parts = split(line);
    parts(parts == "") = [];
    if numel(parts) ~= 13
        continue;
    end
    row = nan(1, 13);
    for j = 1:13
        if parts(j) ~= "-"
            row(j) = str2double(parts(j));
        end
    end
    rows = [rows; row]; %#ok<AGROW>
end

table_data.columns = ["h", "n_L2", "n_L2_order", "n_Linf", "n_Linf_order", ...
    "phi_L2", "phi_L2_order", "phi_Linf", "phi_Linf_order", ...
    "E_L2", "E_L2_order", "E_Linf", "E_Linf_order"];
table_data.data = rows;
end

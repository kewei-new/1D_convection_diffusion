function root = startup_mfem_dd()
%STARTUP_MFEM_DD Add the MFEM-style MATLAB workspace to the path.
root = fileparts(fileparts(mfilename("fullpath")));
addpath(fullfile(root, "matlab"));
addpath(fullfile(root, "cases"));
addpath(fullfile(root, "tests"));
addpath(fullfile(root, "tools"));
end

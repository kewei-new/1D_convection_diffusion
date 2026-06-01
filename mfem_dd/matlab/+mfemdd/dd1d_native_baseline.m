function baseline = dd1d_native_baseline(varargin)
%DD1D_NATIVE_BASELINE Self-contained MATLAB MFEM-DD port of the DD1D MMS run.
opts = local_parse_options(varargin{:});

params = local_default_params();
params.ipdg.method = char(opts.method);
params.dg.p_order = opts.p_order;
params.dg.nloc = opts.p_order + 1;
params.mesh.refine_steps = opts.refine_steps;
params.mesh.ng_base = opts.ng_base;
params.output.verbose = opts.verbose;
params.ipdg = local_method_flags(params.ipdg);
params.ipdg.alpha = (params.dg.p_order + 1)^2;

basis_data = local_basis_data(params);
errors = struct();
errors.n_L2 = zeros(params.mesh.refine_steps, 1);
errors.n_Linf = zeros(params.mesh.refine_steps, 1);
errors.phi_L2 = zeros(params.mesh.refine_steps, 1);
errors.phi_Linf = zeros(params.mesh.refine_steps, 1);
errors.E_L2 = zeros(params.mesh.refine_steps, 1);
errors.E_Linf = zeros(params.mesh.refine_steps, 1);
mesh_sizes = zeros(params.mesh.refine_steps, 1);

for refine_id = 1:params.mesh.refine_steps
    num_cells = params.mesh.ng_base * 2^(refine_id - 1);
    mesh = local_mesh_1d_uniform(params, num_cells);
    mesh_sizes(refine_id) = mesh.cell_sizes(1);
    result = local_solve_1d_dd(mesh, basis_data, params);
    errors.n_L2(refine_id) = local_error_lnorm(result.n_coeff, @local_exact_n, 2, ...
        mesh.num_cells, basis_data.p_order, mesh.cell_sizes, mesh.cell_centers, result.time);
    errors.n_Linf(refine_id) = local_error_absolute(result.n_coeff, @local_exact_n, ...
        mesh.num_cells, basis_data.p_order, mesh.cell_sizes, mesh.cell_centers, result.time);
    errors.phi_L2(refine_id) = local_error_lnorm(result.phi_coeff, @local_exact_phi, 2, ...
        mesh.num_cells, basis_data.p_order, mesh.cell_sizes, mesh.cell_centers, result.time);
    errors.phi_Linf(refine_id) = local_error_absolute(result.phi_coeff, @local_exact_phi, ...
        mesh.num_cells, basis_data.p_order, mesh.cell_sizes, mesh.cell_centers, result.time);
    errors.E_L2(refine_id) = local_error_lnorm(result.E_coeff, @local_exact_phi_x, 2, ...
        mesh.num_cells, basis_data.p_order, mesh.cell_sizes, mesh.cell_centers, result.time);
    errors.E_Linf(refine_id) = local_error_absolute(result.E_coeff, @local_exact_phi_x, ...
        mesh.num_cells, basis_data.p_order, mesh.cell_sizes, mesh.cell_centers, result.time);
end

baseline = struct();
baseline.method = string(params.ipdg.method);
baseline.p_order = params.dg.p_order;
baseline.alpha = params.ipdg.alpha;
baseline.beta = params.ipdg.beta;
baseline.columns = ["h", "n_L2", "n_L2_order", "n_Linf", "n_Linf_order", ...
    "phi_L2", "phi_L2_order", "phi_Linf", "phi_Linf_order", ...
    "E_L2", "E_L2_order", "E_Linf", "E_Linf_order"];
baseline.data = [mesh_sizes, errors.n_L2, local_orders(errors.n_L2), ...
    errors.n_Linf, local_orders(errors.n_Linf), ...
    errors.phi_L2, local_orders(errors.phi_L2), ...
    errors.phi_Linf, local_orders(errors.phi_Linf), ...
    errors.E_L2, local_orders(errors.E_L2), ...
    errors.E_Linf, local_orders(errors.E_Linf)];
baseline.source = "mfem_dd_native_matlab";
baseline.generated_by = "dd1d_native_baseline";
end

function params = local_default_params()
params.domain.left = 0;
params.domain.right = 2 * pi;
params.time.T_end = 1;
params.time.cfl = 0.1;
params.mesh.ng_base = 20;
params.mesh.refine_steps = 5;
params.mesh.uniform = true;
params.dg.p_order = 3;
params.dg.nloc = params.dg.p_order + 1;
params.ipdg.method = "SIPG";
params.ipdg.alpha = (params.dg.p_order + 1)^2;
params.ipdg.beta = 1;
params.model.name = "single_carrier_dd";
params.model.problem_type = "smooth_test";
params.model.dimension = 1;
params.output.verbose = false;
end

function ipdg = local_method_flags(ipdg)
switch upper(string(ipdg.method))
    case "SIPG"
        ipdg.beta = 1;
    case "NIPG"
        ipdg.beta = -1;
    case "IIPG"
        ipdg.beta = 0;
    otherwise
        error("Unknown IPDG method: %s", ipdg.method);
end
end

function basis_data = local_basis_data(params)
basis_data.p_order = params.dg.p_order;
basis_data.nloc = params.dg.nloc;
basis_data.gauss_data = local_gauss_lobatto(params.dg.p_order);
basis_data.inv_mass = local_invmass(params.dg.p_order);
end

function mesh = local_mesh_1d_uniform(params, num_cells)
[cell_sizes, cell_centers] = local_uniform_grid(params.domain.left, ...
    params.domain.right, num_cells);
mesh.left = params.domain.left;
mesh.right = params.domain.right;
mesh.num_cells = num_cells;
mesh.cell_sizes = cell_sizes;
mesh.cell_centers = cell_centers;
mesh.is_uniform = true;
end

function result = local_solve_1d_dd(mesh, basis_data, params)
num_cells = mesh.num_cells;
cell_sizes = mesh.cell_sizes;
cell_centers = mesh.cell_centers;
p_order = basis_data.p_order;
inv_mass = basis_data.inv_mass;
gauss_data = basis_data.gauss_data;
T_end = params.time.T_end;

n_old = local_l2_projection(@local_init_fun, inv_mass, num_cells, ...
    p_order, cell_sizes, cell_centers);

dt = local_evaluate_dt(mesh, params);
diffusion_matrix = local_assemble_ipdg_diffusion_matrix(mesh, basis_data, params);

local_mass_matrix = diag(1 ./ inv_mass) * max(cell_sizes);
mass_matrix = [];
for n = 1:num_cells
    mass_matrix = blkdiag(mass_matrix, local_mass_matrix); %#ok<AGROW>
end

[poisson_matrix_main, poisson_matrix_aux, poisson_rhs_matrix] = ...
    local_assemble_ldg_matrix(gauss_data, inv_mass, p_order, num_cells, cell_sizes);

mass_matrix_dt = mass_matrix / dt;
implicit_matrix = mass_matrix_dt - 0.5 * diffusion_matrix;
time_now = 0;

while time_now < T_end
    if time_now + dt >= T_end
        mass_matrix_dt = mass_matrix_dt * dt;
        dt = T_end - time_now;
        mass_matrix_dt = mass_matrix_dt / dt;
        implicit_matrix = mass_matrix_dt - 0.5 * diffusion_matrix;
    end

    [~, E_coeff] = local_solve_1d_ldg(mesh, basis_data, n_old, time_now, ...
        params, poisson_matrix_main, poisson_matrix_aux, poisson_rhs_matrix);
    transport_rhs_0 = local_transport_rhs(mesh, basis_data, n_old, E_coeff, time_now);
    n_stage1 = implicit_matrix \ (0.5 * transport_rhs_0 + mass_matrix_dt * n_old);

    [~, E_coeff] = local_solve_1d_ldg(mesh, basis_data, n_stage1, ...
        time_now + dt / 2, params, poisson_matrix_main, ...
        poisson_matrix_aux, poisson_rhs_matrix);
    transport_rhs_1 = local_transport_rhs(mesh, basis_data, n_stage1, ...
        E_coeff, time_now + dt / 2);
    n_stage2 = implicit_matrix \ ( ...
        11/18 * transport_rhs_0 + 1/18 * transport_rhs_1 + ...
        1/6 * diffusion_matrix * n_stage1 + mass_matrix_dt * n_old);

    [~, E_coeff] = local_solve_1d_ldg(mesh, basis_data, n_stage2, ...
        time_now + 2 * dt / 3, params, poisson_matrix_main, ...
        poisson_matrix_aux, poisson_rhs_matrix);
    transport_rhs_2 = local_transport_rhs(mesh, basis_data, n_stage2, ...
        E_coeff, time_now + 2 * dt / 3);
    n_stage3 = implicit_matrix \ ( ...
        5/6 * transport_rhs_0 - 5/6 * transport_rhs_1 + 1/2 * transport_rhs_2 + ...
        -1/2 * diffusion_matrix * n_stage1 + 1/2 * diffusion_matrix * n_stage2 + ...
        mass_matrix_dt * n_old);

    [~, E_coeff] = local_solve_1d_ldg(mesh, basis_data, n_stage3, ...
        time_now + dt / 2, params, poisson_matrix_main, ...
        poisson_matrix_aux, poisson_rhs_matrix);
    transport_rhs_3 = local_transport_rhs(mesh, basis_data, n_stage3, ...
        E_coeff, time_now + dt / 2);
    n_new = implicit_matrix \ ( ...
        1/4 * transport_rhs_0 + 7/4 * transport_rhs_1 + ...
        3/4 * transport_rhs_2 - 7/4 * transport_rhs_3 + ...
        3/2 * diffusion_matrix * n_stage1 - 3/2 * diffusion_matrix * n_stage2 + ...
        1/2 * diffusion_matrix * n_stage3 + mass_matrix_dt * n_old);

    n_old = n_new;
    time_now = time_now + dt;
end

[phi_coeff, E_coeff] = local_solve_1d_ldg(mesh, basis_data, n_old, time_now, ...
    params, poisson_matrix_main, poisson_matrix_aux, poisson_rhs_matrix);

result.n_coeff = reshape(n_old, p_order + 1, num_cells);
result.phi_coeff = phi_coeff;
result.E_coeff = E_coeff;
result.time = time_now;
result.dt = dt;
end

function diffusion_matrix = local_assemble_ipdg_diffusion_matrix(mesh, basis_data, params)
p_order = basis_data.p_order;
num_cells = mesh.num_cells;
cell_sizes = mesh.cell_sizes;
gauss_data = basis_data.gauss_data;
alpha = params.ipdg.alpha;
beta = params.ipdg.beta;
Gp = gauss_data(:, 1);
Gw = gauss_data(:, 2);
h = cell_sizes(1);

diffusion_matrix = sparse(num_cells * (p_order + 1), num_cells * (p_order + 1));
Mx1 = zeros(p_order + 1, p_order + 1);
Mx2LL = zeros(p_order + 1, p_order + 1);
Mx2LR = zeros(p_order + 1, p_order + 1);
Mx2RL = zeros(p_order + 1, p_order + 1);
Mx2RR = zeros(p_order + 1, p_order + 1);
Mx3LL = zeros(p_order + 1, p_order + 1);
Mx3LR = zeros(p_order + 1, p_order + 1);
Mx3RR = zeros(p_order + 1, p_order + 1);
Mx3RL = zeros(p_order + 1, p_order + 1);

for d1 = 0:p_order
    for d2 = 0:p_order
        for i1 = 1:length(Gw)
            Mx1(d1 + 1, d2 + 1) = Mx1(d1 + 1, d2 + 1) ...
                + 1 / h * Gw(i1) * local_basis_x(Gp(i1), d2) ...
                * local_basis_x(Gp(i1), d1);
        end
    end
end

for d1 = 0:p_order
    for d2 = 0:p_order
        Mx2LL(d1 + 1, d2 + 1) = 1 / h * local_basis_x(-0.5, d2) * local_basis(-0.5, d1);
        Mx2RR(d1 + 1, d2 + 1) = 1 / h * local_basis_x(0.5, d2) * local_basis(0.5, d1);
        Mx2LR(d1 + 1, d2 + 1) = 1 / h * local_basis_x(-0.5, d2) * local_basis(0.5, d1);
        Mx2RL(d1 + 1, d2 + 1) = 1 / h * local_basis_x(0.5, d2) * local_basis(-0.5, d1);
        Mx3LL(d1 + 1, d2 + 1) = local_basis(-0.5, d2) * local_basis(-0.5, d1);
        Mx3RR(d1 + 1, d2 + 1) = local_basis(0.5, d2) * local_basis(0.5, d1);
        Mx3RL(d1 + 1, d2 + 1) = local_basis(0.5, d2) * local_basis(-0.5, d1);
        Mx3LR(d1 + 1, d2 + 1) = local_basis(-0.5, d2) * local_basis(0.5, d1);
    end
end

Mx2LL = 0.5 * Mx2LL;
Mx2RR = 0.5 * Mx2RR;
Mx2LR = 0.5 * Mx2LR;
Mx2RL = 0.5 * Mx2RL;
Mx3LL = alpha / h * Mx3LL;
Mx3RR = alpha / h * Mx3RR;
Mx3LR = alpha / h * Mx3LR;
Mx3RL = alpha / h * Mx3RL;

index = @(i) (i - 1) * (p_order + 1);
for i = 1:num_cells
    diffusion_matrix(index(i) + 1:index(i + 1), index(i) + 1:index(i + 1)) = ...
        (-Mx1) + ((Mx2RR - Mx2LL) + beta * (Mx2RR' - Mx2LL') - (Mx3RR + Mx3LL));

    if i < num_cells
        diffusion_matrix(index(i) + 1:index(i + 1), index(i + 1) + 1:index(i + 2)) = ...
            (Mx2LR - beta * Mx2RL' + Mx3LR);
    else
        diffusion_matrix(index(1) + 1:index(2), index(i) + 1:index(i + 1)) = ...
            (-Mx2RL + beta * Mx2LR' + Mx3RL);
    end

    if i > 1
        diffusion_matrix(index(i) + 1:index(i + 1), index(i - 1) + 1:index(i)) = ...
            (-Mx2RL + beta * Mx2LR' + Mx3RL);
    else
        diffusion_matrix(index(num_cells) + 1:index(num_cells + 1), index(i) + 1:index(i + 1)) = ...
            (Mx2LR - beta * Mx2RL' + Mx3LR);
    end
end
end

function [M1, M2, Mn] = local_assemble_ldg_matrix(gauss_data, inv_mass, mp, ng, h)
Gp = gauss_data(:, 1);
Gw = gauss_data(:, 2);
M1 = sparse((mp + 1) * ng, (mp + 1) * ng);
C_p = mp / max(h);

for n = 1:ng
    if n == 1
        for alpha = 0:mp
            for beta = 0:mp
                M1(alpha + 1, beta + 1) = inv_mass(alpha + 1) * h(n)^(-1) * ...
                    (-Gw' * (local_basis(Gp, beta) .* local_basis_x(Gp, alpha)) + ...
                    local_basis(1/2, beta) * local_basis(1/2, alpha));
            end
        end
    elseif n < ng
        for alpha = 0:mp
            for beta = 0:mp
                M1((n - 1) * (mp + 1) + (alpha + 1), (n - 1) * (mp + 1) + (beta + 1)) = ...
                    inv_mass(alpha + 1) * h(n)^(-1) * ...
                    (-Gw' * (local_basis(Gp, beta) .* local_basis_x(Gp, alpha)) + ...
                    local_basis(1/2, beta) * local_basis(1/2, alpha));
                M1((n - 1) * (mp + 1) + (alpha + 1), (n - 2) * (mp + 1) + (beta + 1)) = ...
                    inv_mass(alpha + 1) * h(n)^(-1) * ...
                    (-local_basis(1/2, beta) * local_basis(-1/2, alpha));
            end
        end
    else
        for alpha = 0:mp
            for beta = 0:mp
                M1((n - 1) * (mp + 1) + (alpha + 1), (n - 1) * (mp + 1) + (beta + 1)) = ...
                    inv_mass(alpha + 1) * h(n)^(-1) * ...
                    (-Gw' * (local_basis(Gp, beta) .* local_basis_x(Gp, alpha)));
                M1((n - 1) * (mp + 1) + (alpha + 1), (n - 2) * (mp + 1) + (beta + 1)) = ...
                    inv_mass(alpha + 1) * h(n)^(-1) * ...
                    (-local_basis(1/2, beta) * local_basis(-1/2, alpha));
            end
        end
    end
end

M2 = sparse((mp + 1) * ng, (mp + 1) * ng);
M3 = sparse((mp + 1) * ng, (mp + 1) * ng);
for n = 1:ng
    if n < ng
        for alpha = 0:mp
            for beta = 0:mp
                M2((n - 1) * (mp + 1) + (alpha + 1), (n - 1) * (mp + 1) + (beta + 1)) = ...
                    -Gw' * (local_basis(Gp, beta) .* local_basis_x(Gp, alpha)) ...
                    - local_basis(-1/2, beta) * local_basis(-1/2, alpha);
                M2((n - 1) * (mp + 1) + (alpha + 1), n * (mp + 1) + (beta + 1)) = ...
                    local_basis(-1/2, beta) * local_basis(1/2, alpha);
            end
        end
    else
        for alpha = 0:mp
            for beta = 0:mp
                M2((n - 1) * (mp + 1) + (alpha + 1), (n - 1) * (mp + 1) + (beta + 1)) = ...
                    -Gw' * (local_basis(Gp, beta) .* local_basis_x(Gp, alpha)) ...
                    - local_basis(-1/2, beta) * local_basis(-1/2, alpha) ...
                    + local_basis(1/2, beta) * local_basis(1/2, alpha);
                M3((n - 1) * (mp + 1) + (alpha + 1), (n - 1) * (mp + 1) + (beta + 1)) = ...
                    -C_p * local_basis(1/2, beta) * local_basis(1/2, alpha);
            end
        end
    end
end

Mn = M2 * M1 + M3;
end

function [phi_coeff, E_coeff] = local_solve_1d_ldg(mesh, basis_data, carrier_coeff, ...
    time_now, params, M1, M2, Mn)
cell_sizes = mesh.cell_sizes;
cell_centers = mesh.cell_centers;
num_cells = mesh.num_cells;
p_order = basis_data.p_order;
inv_mass = basis_data.inv_mass;
gauss_data = basis_data.gauss_data;
C_p = p_order / max(cell_sizes);
carrier_coeff = reshape(carrier_coeff, p_order + 1, num_cells);
Gp = gauss_data(:, 1);
Gw = gauss_data(:, 2);
f_terms = zeros((p_order + 1) * num_cells, 1);

b1 = zeros(num_cells * (p_order + 1), 1);
n = 1;
for alpha = 0:p_order
    b1(alpha + 1) = -inv_mass(alpha + 1) * cell_sizes(n)^(-1) ...
        * local_ldg_bound(cell_centers(1) - cell_sizes(1) / 2, time_now) ...
        * local_basis(-1/2, alpha);
end

n = num_cells;
for alpha = 0:p_order
    b1((n - 1) * (p_order + 1) + (alpha + 1)) = ...
        inv_mass(alpha + 1) * cell_sizes(n)^(-1) ...
        * local_ldg_bound(cell_centers(num_cells) + cell_sizes(num_cells) / 2, time_now) ...
        * local_basis(1/2, alpha);
end

b2 = zeros(num_cells * (p_order + 1), 1);
for alpha = 0:p_order
    b2((n - 1) * (p_order + 1) + (alpha + 1)) = C_p ...
        * local_ldg_bound(cell_centers(num_cells) + cell_sizes(num_cells) / 2, time_now) ...
        * local_basis(1/2, alpha);
end

for n = 1:num_cells
    Gp_local = Gp * cell_sizes(n) + cell_centers(n);
    f_values = zeros(length(Gp), 1);
    for m = 0:p_order
        f_values = f_values + carrier_coeff(m + 1, n) * local_basis(Gp, m);
    end
    for alpha = 0:p_order
        f_terms((n - 1) * (p_order + 1) + (alpha + 1), 1) = ...
            Gw' * ((f_values - local_f2_fun(Gp_local, time_now)) ...
            .* local_basis(Gp, alpha)) * cell_sizes(n);
    end
end

f_terms = f_terms - M2 * b1 - b2;
phi_coeff_vec = Mn \ f_terms;
E_coeff_vec = M1 * phi_coeff_vec + b1;
phi_coeff = reshape(phi_coeff_vec, p_order + 1, num_cells);
E_coeff = reshape(E_coeff_vec, p_order + 1, num_cells);
end

function transport_rhs = local_transport_rhs(mesh, basis_data, carrier_coeff, E_coeff, time_now)
num_cells = mesh.num_cells;
cell_sizes = mesh.cell_sizes;
cell_centers = mesh.cell_centers;
p_order = basis_data.p_order;
gauss_data = basis_data.gauss_data;
emf = 1;
carrier_coeff = reshape(carrier_coeff, p_order + 1, num_cells);
Gp = gauss_data(:, 1);
Gw = gauss_data(:, 2);
f = zeros(p_order + 1, num_cells);
u1 = ones(length(Gp), num_cells + 2);
phix = ones(length(Gp), num_cells + 2);

for k = 2:num_cells + 1
    u1(:, k) = local_eval_uh(Gp, carrier_coeff(:, k - 1), p_order);
    phix(:, k) = local_eval_uh(Gp, E_coeff(:, k - 1), p_order);
end

u1(2, 1) = u1(2, num_cells + 1);
phix(2, 1) = phix(2, num_cells + 1);
u1(1, num_cells + 2) = u1(1, 2);
phix(1, num_cells + 2) = phix(1, 2);
f1 = u1 .* (-phix);

for k = 2:num_cells + 1
    for m = 0:p_order
        f(m + 1, k - 1) = ...
            -Gw' * (f1(:, k) .* local_basis_x(Gp, m)) ...
            + 0.5 * ((f1(2, k) + f1(1, k + 1)) * local_basis(0.5, m) ...
            - (f1(2, k - 1) + f1(1, k)) * local_basis(-0.5, m)) ...
            - 0.5 * emf * ((u1(1, k + 1) - u1(2, k)) * local_basis(0.5, m) ...
            - (u1(1, k) - u1(2, k - 1)) * local_basis(-0.5, m));
    end
end

hu = zeros(p_order + 1, num_cells);
for k = 2:num_cells + 1
    for m = 0:p_order
        hu(m + 1, k - 1) = ...
            Gw' * (local_f1_fun(Gp * cell_sizes(k - 1) + cell_centers(k - 1), time_now) ...
            .* local_basis(Gp, m)) * cell_sizes(1);
    end
end

transport_rhs = reshape(f, num_cells * (p_order + 1), 1) + ...
    reshape(hu, num_cells * (p_order + 1), 1);
end

function uh0 = local_l2_projection(initial_fun, inv_a, ng, mp, h, mid_points)
gauss_data = local_gauss_legendre(4);
Gp = gauss_data(:, 1);
Gw = gauss_data(:, 2);
uh0 = zeros(mp + 1, ng);
for n = 1:ng
    Gp_local = Gp * h(n) + mid_points(n);
    ini_value = initial_fun(Gp_local);
    for m = 0:mp
        uh0(m + 1, n) = inv_a(m + 1) * Gw' * (ini_value .* local_basis(Gp, m));
    end
end
uh0 = reshape(uh0, (mp + 1) * ng, 1);
end

function e = local_error_lnorm(uh, exact_fun, L_order, ng, mp, h, mid_points, t)
if isvector(uh)
    uh = reshape(uh, mp + 1, ng);
elseif size(uh, 1) ~= mp + 1
    uh = uh.';
end
gauss_data = local_gauss_lobatto(mp);
Gp = gauss_data(:, 1);
Gw = gauss_data(:, 2);
e = 0;
for n = 1:ng
    x = Gp * h(n) + mid_points(n);
    uh_value = local_eval_uh(Gp, uh(:, n), mp);
    uex = exact_fun(x, t);
    if L_order == 2
        e = e + sum(Gw .* (uh_value - uex).^2) * h(n);
    else
        e = e + sum(Gw .* abs(uh_value - uex).^L_order) * h(n);
    end
end
e = e^(1 / L_order);
end

function e = local_error_absolute(uh, exact_fun, ng, mp, h, mid_points, t)
if isvector(uh)
    uh = reshape(uh, mp + 1, ng);
elseif size(uh, 1) ~= mp + 1
    uh = uh.';
end
gauss_data = local_gauss_lobatto(mp);
Gp = gauss_data(:, 1);
e = 0;
for n = 1:ng
    x = Gp * h(n) + mid_points(n);
    uh_value = local_eval_uh(Gp, uh(:, n), mp);
    uex = exact_fun(x, t);
    e = max(e, max(abs(uh_value - uex)));
end
end

function result = local_eval_uh(x, uh, mp)
result = 0;
for m = 0:mp
    result = result + uh(m + 1) * local_basis(x, m);
end
end

function dt = local_evaluate_dt(mesh, params)
dt = params.time.cfl * mesh.cell_sizes(1);
end

function [h, mid_points] = local_uniform_grid(left, right, ng)
h0 = (right - left) / ng;
mid_points = linspace(left + h0 / 2, right - h0 / 2, ng);
h = ones(1, ng) * h0;
end

function inv_a = local_invmass(mp)
if mp == 1
    inv_a = [1, 12];
elseif mp == 2
    inv_a = [1, 12, 180];
elseif mp == 3
    inv_a = [1, 12, 180, 2800];
else
    error("Unsupported polynomial order for DD1D native backend: %d", mp);
end
end

function gauss = local_gauss_lobatto(mp)
if mp == 1
    gauss = [-0.5, 1.0/6.0; 0.5, 1.0/6.0; 0.0, 2.0/3.0];
elseif mp == 2
    gauss = [-0.5, 1.0/12.0; 0.5, 1.0/12.0; ...
        -sqrt(5.0)/10.0, 5.0/12.0; sqrt(5.0)/10.0, 5.0/12.0];
elseif mp == 3
    gauss = [-0.5, 1.0/20.0; 0.5, 1.0/20.0; ...
        -sqrt(21.0)/14.0, 49.0/180.0; ...
        sqrt(21.0)/14.0, 49.0/180.0; 0.0, 64.0/180.0];
elseif mp == 4
    gauss = zeros(6, 2);
    gauss(1, :) = [-0.5, 1.0/30.0];
    gauss(2, :) = [0.5, 1.0/30.0];
    gauss(3, 1) = -sqrt(147.0 + 42.0 * sqrt(7.0)) / 42.0;
    gauss(4, 1) = sqrt(147.0 + 42.0 * sqrt(7.0)) / 42.0;
    gauss(5, 1) = -sqrt(147.0 - 42.0 * sqrt(7.0)) / 42.0;
    gauss(6, 1) = sqrt(147.0 - 42.0 * sqrt(7.0)) / 42.0;
    gauss(3, 2) = (-7.0 + 5.0 * sqrt(7.0)) * sqrt(7.0) * (7.0 + sqrt(7.0)) / 840.0;
    gauss(4, 2) = gauss(3, 2);
    gauss(5, 2) = (7.0 + 5.0 * sqrt(7.0)) * sqrt(7.0) / (7.0 + sqrt(7.0)) / 20.0;
    gauss(6, 2) = gauss(5, 2);
else
    error("Unsupported Gauss-Lobatto order: %d", mp);
end
end

function gauss = local_gauss_legendre(n)
switch n
    case 4
        x_std = [-0.8611363115940526; -0.3399810435848563; ...
            0.3399810435848563; 0.8611363115940526];
        w_std = [0.3478548451374538; 0.6521451548625461; ...
            0.6521451548625461; 0.3478548451374538];
    otherwise
        error("Unsupported Gauss-Legendre order: %d", n);
end
gauss = [x_std / 2, w_std / 2];
end

function res = local_basis(x, mp)
if mp == 0
    res = 0 .* x + 1;
elseif mp == 1
    res = 1 .* x;
elseif mp == 2
    res = x.^2 - 1 ./ 12;
elseif mp == 3
    res = x.^3 - 0.15 * x;
elseif mp == 4
    res = (x.^2 - 3 / 14) .* x .* x + 3.0 / 560;
else
    error("Unsupported basis order: %d", mp);
end
end

function res = local_basis_x(x, mp)
if mp == 0
    res = 0 .* x;
elseif mp == 1
    res = 0 .* x + 1;
elseif mp == 2
    res = 2 .* x;
elseif mp == 3
    res = 3 * x.^2 - 0.15;
elseif mp == 4
    res = 4 * x.^3 - 3 / 7 * x;
else
    error("Unsupported basis derivative order: %d", mp);
end
end

function result = local_exact_n(x, t)
result = sin(x) * cos(t);
end

function result = local_exact_phi(x, t)
result = sin(t) * cos(x);
end

function result = local_exact_phi_x(x, t)
result = -sin(t) * sin(x);
end

function result = local_f1_fun(x, t)
result = -sin(t) .* sin(x) + cos(t) .* sin(x) ...
    - 2 * sin(t) .* cos(t) .* sin(x) .* cos(x);
end

function result = local_f2_fun(x, t)
result = sin(t) * cos(x) + sin(x) * cos(t);
end

function result = local_init_fun(x)
result = sin(x);
end

function result = local_ldg_bound(x, t)
if abs(x) < 1.0e-10 || abs(x - 2 * pi) < 1.0e-10
    result = sin(t);
else
    error("LDG boundary requested away from a DD1D endpoint: %.16g", x);
end
end

function orders = local_orders(errors)
orders = nan(size(errors));
for k = 2:numel(errors)
    if errors(k) > 0 && errors(k - 1) > 0
        orders(k) = log2(errors(k - 1) / errors(k));
    end
end
end

function opts = local_parse_options(varargin)
opts = struct("method", "SIPG", "p_order", 3, ...
    "refine_steps", 5, "ng_base", 20, "verbose", false);
if mod(numel(varargin), 2) ~= 0
    error("Options must be name/value pairs.");
end
for k = 1:2:numel(varargin)
    name = lower(string(varargin{k}));
    value = varargin{k + 1};
    switch name
        case "method"
            opts.method = string(value);
        case "p_order"
            opts.p_order = value;
        case "refine_steps"
            opts.refine_steps = value;
        case "ng_base"
            opts.ng_base = value;
        case "verbose"
            opts.verbose = logical(value);
        otherwise
            error("Unknown option: %s", name);
    end
end
end

function mesh = build_mesh_1D_uniform(params, num_cells)
%BUILD_MESH_1D_UNIFORM Build uniform 1D mesh struct.

    [cell_sizes, cell_centers] = generate_1D_uniform_grid( ...
        params.domain.left, params.domain.right, num_cells);

    mesh.left         = params.domain.left;
    mesh.right        = params.domain.right;
    mesh.num_cells    = num_cells;
    mesh.cell_sizes   = cell_sizes;
    mesh.cell_centers = cell_centers;
    mesh.is_uniform   = true;
end
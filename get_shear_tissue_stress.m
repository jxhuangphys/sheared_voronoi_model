function mean_stress_int = get_shear_tissue_stress(center_xy, cell_chain, vertex_position, ...
                K_A, A0_list, K_P, P0_list, gamma, box_size)
% This function calculates the interaction stress tensor in a sheared tissue.
%
% Inputs:
% - center_xy: Positions of cell centers [N_cell x 2].
% - cell_chain: Cell adjacency information, list of vertex indices for each cell.
% - vertex_position: Coordinates of all vertices in the Voronoi tessellation.
% - K_A: Area elasticity constant for each cell [N_cell x 1].
% - A0_list: Preferred area for each cell [N_cell x 1].
% - K_P: Perimeter elasticity constant for each cell [N_cell x 1].
% - P0_list: Preferred perimeter for each cell [N_cell x 1].
% - gamma: Shear strain applied to the system.
% - box_size: Size of the simulation box.
%
% Outputs:
% - mean_stress_int: Mean interaction stress tensor [1 x 4] for the system.

% Normalize gamma to prevent overflow in cyclic boundaries
gamma = mod(gamma, 1);

% Initialize key parameters
N_cell = numel(cell_chain); % Number of cells in the system
if ~exist('box_size', 'var') % If box size not provided, calculate from the number of cells
    box_size = sqrt(N_cell);
end
Lx = box_size(1); % Width of the simulation box
Ly = box_size(end); % Height of the simulation box

% Initialize lists to store perimeter and area for each cell
perimeter_list = zeros(N_cell, 1);
area_list = zeros(N_cell, 1);

% Loop over each cell to calculate perimeter and area
for i_c = 1:N_cell
    chain_tail2 = cell_chain{i_c}; % Vertices forming the current cell
    chain_tail2 = [chain_tail2, chain_tail2(1:2)]; % Close the loop for periodic boundaries
    xy_tail2 = vertex_position(chain_tail2, :); % Vertex positions for the current cell

    % Apply shear corrections for periodic boundary conditions
    if center_xy(i_c, 2) > Ly - 3
        xy_tail2(xy_tail2(:, 2) <= 3, 1) = xy_tail2(xy_tail2(:, 2) <= 3, 1) + gamma * Lx;
    elseif center_xy(i_c, 2) < 3
        xy_tail2(xy_tail2(:, 2) >= Ly - 3, 1) = xy_tail2(xy_tail2(:, 2) >= Ly - 3, 1) - gamma * Lx;
    end

    % Relocate vertices to account for periodic boundary conditions
    xy_tail2 = pbc_relocate(center_xy(i_c, :), xy_tail2, box_size);

    % Compute the perimeter by summing the edge lengths
    dx_dy_tail1 = xy_tail2(2:end-1, :) - xy_tail2(1:end-2, :);
    perimeter_list(i_c) = sum(vecnorm(dx_dy_tail1, 2, 2));

    % Compute the area using the shoelace formula
    dy2 = xy_tail2(3:end, 2) - xy_tail2(1:end-2, 2);
    area_list(i_c) = abs(sum(xy_tail2(2:end-1, 1) .* dy2)) / 2;
end

% Obtain adjacency information for cells and edges
[~, cell_edge_adj, edge_vid] = get_tissue_adj(cell_chain);

% Initialize edge vectors
N_edge = size(edge_vid, 1); % Number of edges
edge_vector = zeros(N_edge, 2); % Edge vectors initialization

% Loop over edges to calculate edge vectors
for i_e = 1:N_edge
    two_vid = edge_vid(i_e, :); % Indices of the two vertices forming the edge
    two_end_position = vertex_position(two_vid, :); % Positions of the two vertices

    % Correct for periodic boundary conditions in the y-direction
    if two_end_position(1, 2) - two_end_position(2, 2) > Ly / 2
        two_end_position(1, 1) = two_end_position(1, 1) - Lx * gamma;
    elseif two_end_position(2, 2) - two_end_position(1, 2) > Ly / 2
        two_end_position(1, 1) = two_end_position(1, 1) + Lx * gamma;
    end

    % Relocate vertices for periodic boundary conditions
    two_end_position = pbc_relocate(two_end_position(1, :), two_end_position, box_size);

    % Calculate the edge vector
    edge_vector(i_e, :) = two_end_position(1, :) - two_end_position(2, :);
end

% Compute edge lengths and initialize tension
edge_len = vecnorm(edge_vector, 2, 2);

% Compute cell perimeter deviations and edge tensions
cell_perimeter_list = cell_edge_adj * edge_len;
cell_peri_devia_list = 2 * K_P .* (cell_perimeter_list - P0_list);
edge_tension_list = cell_edge_adj.' * cell_peri_devia_list;

% Compute cell pressures based on area deviations
cell_pressure_list = -2 * K_A .* (area_list - A0_list);

% Eliminate numerical issues for very short edges
short_edge_id = find(edge_len < 1e-15);
edge_len(short_edge_id) = 1; % Avoid division by zero
edge_vector(short_edge_id, :) = 0; % Zero out short edges

% Initialize interaction stress tensor
stress_int = zeros(N_cell, 4); % Tensor components [xx, xy, yx, yy]

% Loop over cells to calculate stress tensor components
for i_c = 1:N_cell
    i_edge = find(cell_edge_adj(i_c, :)); % Edges associated with the current cell
    for alpha = 1:2
        for beta = alpha:2
            col_ind = (alpha - 1) * 2 + beta; % Column index for stress tensor
            stress_int(i_c, col_ind) = -cell_pressure_list(i_c) * (alpha == beta) ...
                + 1 / (2 * area_list(i_c)) * sum(edge_tension_list(i_edge) .* ...
                edge_vector(i_edge, alpha) .* edge_vector(i_edge, beta) ./ edge_len(i_edge));
        end
    end
end

% Symmetrize the stress tensor components
stress_int(:, 3) = stress_int(:, 2);

% Calculate mean interaction stress tensor for the system
mean_stress_int = sum(area_list .* stress_int, 1) / N_cell;

end % End of the function

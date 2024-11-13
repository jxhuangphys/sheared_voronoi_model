function [center_force, dE_dH_mat] = ...
    get_shear_voronoi_force(center_xy, cell_chain, vertex_position, ...
                            gamma, K_A, A0_list, K_P, P0_list, area_list, perimeter_list)

% Number of cells and vertices in the tissue
N_cell = size(center_xy, 1); % Number of cells
N_vertex = size(vertex_position, 1); % Number of vertices in the Voronoi tessellation
box_size = sqrt(N_cell); % Calculate box size based on cell number

% Apply Lees–Edwards boundary condition to gamma
gamma = mod(gamma, 1); % Gamma is reduced to the interval [0, 1)
if gamma > 0.5
    gamma = gamma - 1; % Shift gamma to the interval [-0.5, 0.5] for symmetry
end

% Initialize adjacency matrix and energy gradient storage
adj_cell_vertex = zeros(N_cell, N_vertex); % Tracks cell-vertex adjacency
dE_dH_mat = zeros(N_cell, 2 * N_vertex); % Gradient of energy with respect to vertex positions

% Loop over all cells to compute energy gradients
for i_c = 1:N_cell
    chain = cell_chain{i_c}; % Get vertices of the current cell
    chain_len = numel(chain); % Number of vertices in the cell
    adj_cell_vertex(i_c, chain) = 1; % Update adjacency matrix
    chain = [chain, chain(1), chain(2)]; % Make vertex list cyclic for derivative calculation

    % Get vertex positions (x, y) and adjust for shear-induced displacements
    x_list = vertex_position(chain, 1); 
    y_list = vertex_position(chain, 2);

    % Apply Lees–Edwards boundary condition for shear deformation
    if center_xy(i_c, 2) > box_size * 2 / 3
        % Handle cells near the upper edge
        x_list(y_list <= box_size / 3) = x_list(y_list <= box_size / 3) + gamma * box_size;
    elseif center_xy(i_c, 2) < box_size / 3
        % Handle cells near the lower edge
        x_list(y_list >= box_size * 2 / 3) = x_list(y_list >= box_size * 2 / 3) - gamma * box_size;
    end

    % Apply periodic boundary conditions (PBC) to ensure all vertices are within the simulation box
    x_list = pbc_relocate(center_xy(i_c, 1), x_list, box_size);
    y_list = pbc_relocate(center_xy(i_c, 2), y_list, box_size);

    % Compute differences between neighboring vertices for energy calculations
    dx_list = diff(x_list); % Differences in x-coordinates (v+1 edges)
    dy_list = diff(y_list); % Differences in y-coordinates (v+1 edges)
    norm_list = sqrt(dx_list.^2 + dy_list.^2); % Norm of edge vectors
    d2x_list = x_list(3:end) - x_list(1:end-2); % Second differences in x (v edges)
    d2y_list = y_list(3:end) - y_list(1:end-2); % Second differences in y (v edges)

    % Handle overlapping vertices by adjusting zero-length edges
    if ~all(norm_list)
        index_list = find(norm_list == 0); % Identify overlapping edges
        norm_list(index_list) = 1; % Avoid division by zero
        for i_index = 1:numel(index_list)
            id = index_list(i_index);
            left_id = mod(id - 2, chain_len) + 1;
            right_id = mod(id, chain_len) + 1;
            if dx_list(left_id) * dx_list(right_id) > 0
                dx_list(id) = sign(dx_list(left_id));
            end
            if dy_list(left_id) * dy_list(right_id) > 0
                dy_list(id) = sign(dy_list(left_id));
            end
        end
    end

    % Calculate energy gradients with respect to vertex positions
    for n = 1:chain_len
        % Gradient w.r.t. x-coordinates of vertices
        dE_dH_mat(i_c, chain(n+1)) = ...
            K_A(i_c) * (area_list(i_c) - A0_list(i_c)) * d2y_list(n) + ...
            2 * K_P(i_c) * (perimeter_list(i_c) - P0_list(i_c)) * dx_list(n) / norm_list(n) - ...
            2 * K_P(i_c) * (perimeter_list(i_c) - P0_list(i_c)) * dx_list(n+1) / norm_list(n+1);

        % Gradient w.r.t. y-coordinates of vertices
        dE_dH_mat(i_c, chain(n+1) + N_vertex) = ...
            K_A(i_c) * (area_list(i_c) - A0_list(i_c)) * (-d2x_list(n)) + ...
            2 * K_P(i_c) * (perimeter_list(i_c) - P0_list(i_c)) * dy_list(n) / norm_list(n) - ...
            2 * K_P(i_c) * (perimeter_list(i_c) - P0_list(i_c)) * dy_list(n+1) / norm_list(n+1);
    end
end

% Calculate partial derivatives of vertex positions w.r.t. cell center positions
dh_dc_mat = zeros(2 * N_vertex, 2 * N_cell);
% calculations for dh_dc_mat (defined in Bi 2016 PRX Appendix A for force calculations)
for i_v = 1:N_vertex

	joint_v_x = vertex_position(i_v, 1);
	joint_v_y = vertex_position(i_v, 2);

	tri_vid = find(adj_cell_vertex(:, i_v)>0);
	e1_id = tri_vid(1);
	e2_id = tri_vid(2);
	e3_id = tri_vid(3);
    
    triple_center_xy = center_xy(tri_vid, :);
    % undo lee edwards shift
    if joint_v_y>box_size*2/3
        triple_center_xy(triple_center_xy(:, 2)<=box_size/3, 1) = ...
            triple_center_xy(triple_center_xy(:, 2)<=box_size/3, 1) ...
            + gamma*box_size;
    elseif joint_v_y<box_size/3
        triple_center_xy(triple_center_xy(:, 2)>=box_size*2/3, 1) = ...
            triple_center_xy(triple_center_xy(:, 2)>=box_size*2/3, 1) ...
            - gamma*box_size;
    end

    x1 = pbc_relocate(joint_v_x, triple_center_xy(1, 1), box_size);
    y1 = pbc_relocate(joint_v_y, triple_center_xy(1, 2), box_size);
    x2 = pbc_relocate(joint_v_x, triple_center_xy(2, 1), box_size);
    y2 = pbc_relocate(joint_v_y, triple_center_xy(2, 2), box_size);
    x3 = pbc_relocate(joint_v_x, triple_center_xy(3, 1), box_size);
    y3 = pbc_relocate(joint_v_y, triple_center_xy(3, 2), box_size);

    % d hx / d e1x
    dh_dc_mat(i_v, e1_id) = -((y2-y3)*(2*x1*(x3*(y2-y1)+x2*(y1-y3))+(y1-y2)*(x3^2+(y1-y3)*(y2-y3))+x2^2*(y3-y1)+ x1^2*(y3-y2)))/(2*(x3*(y2-y1)+x2*(y1-y3)+x1*(y3-y2))^2);
    % d hx / d e2x
    dh_dc_mat(i_v, e2_id) =  ((x3^2*(y1-y2)-2*x2*(x3*(y1-y2)+x1*(y2-y3))+x2^2*(y1-y3)+(x1^2+(y1-y2)*(y1-y3))*(y2-y3))*(y1-y3))/(2*(x3*(y2-y1)+x2*(y1-y3)+x1*(y3-y2))^2);
    % d hx / d e3x
    dh_dc_mat(i_v, e3_id) =  -((y1-y2)*(x3^2*(y2-y1)+2*x2*x3*(y1-y3)+(x1^2+(y1-y2)*(y1-y3))*(y2-y3)+x2^2*(y3-y1)+ 2*x1*x3*(y3-y2)))/(2*(x3*(y2-y1)+x2*(y1-y3)+x1*(y3-y2))^2);

    % d hx / d e1y
    dh_dc_mat(i_v, e1_id + N_cell) = ((y2-y3)*(x1^2*(x2-x3)+x2^2*x3+x3*(y1-y2)^2-x2*(x3^2+(y1-y3)^2)+ x1*(-x2^2+x3^2+2*y1*y2-y2^2-2*y1*y3+y3^2)))/(2*(x3*(y2-y1)+x2*(y1-y3)+x1*(y3-y2))^2);
    % d hx / d e2y
    dh_dc_mat(i_v, e2_id + N_cell) = -((y1-y3)*(x1^2*(x2-x3)+x2^2*x3-x2*x3^2-x3*(y1-y2)^2+x1*(-x2^2+x3^2+(y2-y3)^2)+ x2*(y1-y3)*(y1-2*y2+y3)))/(2*(x3*(y2-y1)+x2*(y1-y3)+x1*(y3-y2))^2);
    % d hx / d e3y
    dh_dc_mat(i_v, e3_id + N_cell) = ((y1-y2)*(x1^2*(x2-x3)+x2^2*x3+x2*(-x3^2+(y1-y3)^2)-x1*(x2^2-x3^2+(y2-y3)^2)- x3*(y1-y2)*(y1+y2-2*y3)))/(2*(x3*(y2-y1)+x2*(y1-y3)+x1*(y3-y2))^2);

    % d hy / d e1x
    dh_dc_mat(i_v + N_vertex, e1_id) = -((x2-x3)*(-((y1-y2)*(x3^2+(y1-y3)*(y2-y3)))+x2^2*(y1-y3)+x1^2*(y2-y3)+ 2*x1*(x3*(y1-y2)+x2*(y3-y1))))/(2*(x3*(y2-y1)+x2*(y1-y3)+x1*(y3-y2))^2);
    % d hy / d e2x
    dh_dc_mat(i_v + N_vertex, e2_id) = -((x1-x3)*(x3^2*(y1-y2)-2*x2*(x3*(y1-y2)+x1*(y2-y3))+x2^2*(y1-y3)+ (x1^2+(y1-y2)*(y1-y3))*(y2-y3)))/(2*(x3*(y2-y1)+x2*(y1-y3)+x1*(y3-y2))^2);
    % d hy / d e3x
    dh_dc_mat(i_v + N_vertex, e3_id) = ((x1-x2)*(x3^2*(y2-y1)+2*x2*x3*(y1-y3)+(x1^2+(y1-y2)*(y1-y3))*(y2-y3)+x2^2*(y3-y1)+ 2*x1*x3*(y3-y2)))/(2*(x3*(y2-y1)+x2*(y1-y3)+x1*(y3-y2))^2);

    % d hy / d e1y
    dh_dc_mat(i_v + N_vertex, e1_id + N_cell) = -((x2-x3)*(x1^2*(x2-x3)+x2^2*x3+x3*(y1-y2)^2-x2*(x3^2+(y1-y3)^2)+ x1*(-x2^2+x3^2+2*y1*y2-y2^2-2*y1*y3+y3^2)))/(2*(x3*(y2-y1)+x2*(y1-y3)+x1*(y3-y2))^2);
    % d hy / d e2y
    dh_dc_mat(i_v + N_vertex, e2_id + N_cell) = ((x1-x3)*(x1^2*(x2-x3)+x2^2*x3-x2*x3^2-x3*(y1-y2)^2+x1*(-x2^2+x3^2+(y2-y3)^2)+ x2*(y1-y3)*(y1-2*y2+y3)))/(2*(x3*(y2-y1)+x2*(y1-y3)+x1*(y3-y2))^2);
    % d hy / d e3y
    dh_dc_mat(i_v + N_vertex, e3_id + N_cell) = -((x1-x2)*(x1^2*(x2-x3)+x2^2*x3+x2*(-x3^2+(y1-y3)^2)-x1*(x2^2-x3^2+(y2-y3)^2)- x3*(y1-y2)*(y1+y2-2*y3)))/(2*(x3*(y2-y1)+x2*(y1-y3)+x1*(y3-y2))^2);

	
end


% Compute forces on cell centers by combining energy gradients and derivatives
center_force = zeros(N_cell, 2);
dH_mat_x = dh_dc_mat(:, 1:N_cell); % x-component of derivative matrix
dH_mat_y = dh_dc_mat(:, N_cell + (1:N_cell)); % y-component of derivative matrix
center_force(:, 1) = -sum(dE_dH_mat * dH_mat_x, 1);
center_force(:, 2) = -sum(dE_dH_mat * dH_mat_y, 1);

end

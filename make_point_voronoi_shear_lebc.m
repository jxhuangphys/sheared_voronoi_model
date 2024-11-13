function [cell_chain, vertex_position, area_list, perimeter_list] = ...
    make_point_voronoi_shear_lebc(center_xy, gamma)

% Function to compute Voronoi tessellation under shear deformation
% with Lees-Edwards boundary conditions (LEBC).
% This function returns the Voronoi cell chains, vertex positions,
% cell areas, and perimeters.
%
% Inputs:
% - center_xy: Positions of cell centers
% - gamma: Shear strain applied
%
% Outputs:
% - cell_chain: List of vertex indices defining each cell
% - vertex_position: Coordinates of the vertices in the central box
% - area_list: List of cell areas
% - perimeter_list: List of cell perimeters

% Number of cells and simulation box size
N_cell = size(center_xy, 1);
box_size = sqrt(N_cell);

% Apply periodicity to the shear strain
gamma = mod(gamma, 1);
if gamma > 0.5 
    gamma = gamma - 1;
end

% Shift positions into the central simulation box
center_xy = mod(center_xy, box_size);

% Extend the simulation to include neighboring layers for proper tessellation
extend_layer = min(5, 3 + sqrt(N_cell / 100)); % At least 3 layers, at most 5 layers

% Duplicate positions for Lees-Edwards boundary conditions (LEBC)
% Add duplicates for left/right and top/bottom with shear applied
center_xy = [center_xy; center_xy + box_size * [1 0]; ...
             center_xy - box_size * [1 0]]; % Extend horizontally
center_xy = [center_xy; center_xy + box_size * [gamma 1]; ...
             center_xy - box_size * [gamma 1]]; % Extend vertically with shear

% Restrict positions to within the extended simulation boundary
center_xy = center_xy(all(center_xy > -extend_layer, 2) & ...
                      all(center_xy < box_size + extend_layer, 2), :);

% Generate Voronoi diagram using Delaunay triangulation
DT = delaunayTriangulation(center_xy);
[big_vertex_position, big_cell_chain] = voronoiDiagram(DT);

% Extract cells within the central box
cell_chain = big_cell_chain(1:N_cell);
central_to_big_vid_dict = find(all(big_vertex_position >= 0, 2) & ...
                               all(big_vertex_position < box_size, 2));

% Create mapping from global to central vertices
big_to_central_vid_dict = zeros(1, size(big_vertex_position, 1));
big_to_central_vid_dict(central_to_big_vid_dict) = 1:numel(central_to_big_vid_dict);

% Restrict vertices to those inside the central box
vertex_position = mod(big_vertex_position(central_to_big_vid_dict, :), box_size);

% Handle vertices outside the central box due to shear
extend_pool = unique([cell_chain{:}]);
outside_vid = setdiff(extend_pool, central_to_big_vid_dict);
shifted_ov_xy = big_vertex_position(outside_vid, :);

% Correct horizontal positions for shear-induced displacements
shifted_ov_xy(:, 1) = shifted_ov_xy(:, 1) - box_size * gamma * ...
    ((shifted_ov_xy(:, 2) > box_size) - (shifted_ov_xy(:, 2) < 0));
shifted_ov_xy = mod(shifted_ov_xy, box_size);

% Match outside vertices to nearest central vertices
for i_o = 1:numel(outside_vid)
    i_o_xy = shifted_ov_xy(i_o, :);
    [~, min_pos] = min(sum(abs(i_o_xy - vertex_position), 2));
    big_to_central_vid_dict(outside_vid(i_o)) = min_pos;
end

% Map cell chains to use central vertex indices
cell_chain = cellfun(@(c) big_to_central_vid_dict(c), cell_chain, 'uniformoutput', 0);

% Check for rosette defects (vertices with incorrect connectivity)
if min(diff(sort(sum(vertex_position, 2)))) < 1e-13
    % Handle vertices with incorrect connectivity
    cell_vertex_adj = zeros(N_cell, size(vertex_position, 1));
    for i_c = 1:N_cell
        cell_vertex_adj(i_c, cell_chain{i_c}) = 1;
    end
    
    % Correct cells with missing or duplicate connections
    lack_connect_vid = find(sum(cell_vertex_adj, 1) < 3, 1);
    if ~isempty(lack_connect_vid)
        cell_repeat_ind = find(cellfun(@(c) numel(c) - numel(unique(c)), cell_chain));
        % Correct overlapping vertices within cell chains
        for i_wc = 1:numel(cell_repeat_ind)
            wrong_cell = cell_repeat_ind(i_wc);
            chain = cell_chain{wrong_cell};
            v_small_id = mode(chain);
            v_big_id = find(sum(abs(vertex_position - vertex_position(v_small_id, :)), 2) < 1e-13);
            v_big_id(v_big_id == v_small_id) = [];
            v_big_id = v_big_id(ismember(v_big_id, lack_connect_vid));
            chain(find(chain == v_small_id, 1, 'first')) = v_big_id(1);
            cell_chain{wrong_cell} = chain;
        end
    end
end

% Initialize area and perimeter lists
perimeter_list = zeros(N_cell, 1);
area_list = zeros(N_cell, 1);

% Compute area and perimeter for each cell
for i_c = 1:N_cell
    chain_xy = vertex_position(cell_chain{i_c}, :);
    chain_len = size(chain_xy, 1);

    % Handle shear-induced shifts near the box edges
    if center_xy(i_c, 2) > box_size - 3
        chain_xy(chain_xy(:, 2) <= 3, 1) = chain_xy(chain_xy(:, 2) <= 3, 1) + gamma * box_size;
    elseif center_xy(i_c, 2) < 3
        chain_xy(chain_xy(:, 2) >= box_size - 3, 1) = chain_xy(chain_xy(:, 2) >= box_size - 3, 1) - gamma * box_size;
    end

    % Relocate vertices relative to the cell center
    chain_xy = pbc_relocate(center_xy(i_c, :), chain_xy, box_size);

    % Compute edge vectors
    dx = chain_xy([2:chain_len 1], 1) - chain_xy(:, 1);
    dy = chain_xy([2:chain_len 1], 2) - chain_xy(:, 2);

    % Calculate area using boundary integral
    area_list(i_c) = -sum(chain_xy(:, 2) .* dx - chain_xy(:, 1) .* dy) / 2;

    % Calculate perimeter as the sum of edge lengths
    perimeter_list(i_c) = sum(sqrt(dx.^2 + dy.^2));
end

end

function [cell_cell_adj, cell_edge_adj, edge_vid, ...
          edge_vector, edge_tension_list] = ...
          get_tissue_adj(cell_chain, vertex_position, box_size, K_P, P0_list)
% Function to compute various adjacency and geometric relationships
% within a Voronoi-based tissue model.
% Inputs:
%   - cell_chain: Cell arrays containing vertex indices for each cell.
%   - vertex_position: Coordinates of vertices in the simulation box.
%   - box_size: Size of the periodic simulation box.
%   - K_P: Perimeter elasticity constant for cells.
%   - P0_list: List of preferred perimeters for cells.
% Outputs:
%   - cell_cell_adj: Adjacency matrix indicating neighboring cells.
%   - cell_edge_adj: Matrix linking cells to shared edges.
%   - edge_vid: Matrix with vertex IDs defining each edge.
%   - edge_vector: Vector representing edge lengths and directions.
%   - edge_tension_list: List of tensions for each edge.

N_cell = numel(cell_chain); % Total number of cells
if nargin == 1 && nargout > 3
    error('Vertex positions required for edge length calculations');
end

% Generate the cell-vertex adjacency matrix
cell_vertex_adj = zeros(N_cell, 2 * N_cell); % Preallocate adjacency matrix
for i_c = 1:N_cell
    cell_vertex_adj(i_c, cell_chain{i_c}) = 1; % Fill adjacency for each cell
end

% Compute the cell-cell adjacency matrix
cell_cell_adj = zeros(N_cell, N_cell); % Preallocate adjacency matrix
for i_v = 1:2 * N_cell % Loop over all potential vertices
    triple_cell = find(cell_vertex_adj(:, i_v) > 0); % Cells sharing a vertex
    if isempty(triple_cell), continue, end % Skip if no cells share this vertex
    if numel(triple_cell) > 1 % If two cells share the vertex
        cell_cell_adj(triple_cell(1), triple_cell(2)) = 1;
    end
    if numel(triple_cell) > 2 % If three cells share the vertex
        cell_cell_adj(triple_cell(2), triple_cell(3)) = 1;
        cell_cell_adj(triple_cell(3), triple_cell(1)) = 1;
    end
end
cell_cell_adj = (cell_cell_adj + cell_cell_adj.') > 0; % Symmetrize and convert to boolean

% If additional outputs are required, compute cell-edge adjacency
if nargout >= 2
    % Find all unique edges based on cell adjacency
    ids = find(cell_cell_adj > 0); % Find adjacent cell pairs
    [a, b] = ind2sub(size(cell_cell_adj), ids); % Get indices of the pairs
    edgelist = unique(sort([a, b], 2), 'rows'); % Create unique edge list
    N_edge = size(edgelist, 1); % Total number of edges

    % Create cell-edge adjacency matrix
    cell_edge_adj = zeros(N_cell, N_edge); % Preallocate matrix
    cell_edge_adj(sub2ind(size(cell_edge_adj), edgelist(:, 1)', 1:N_edge)) = 1;
    cell_edge_adj(sub2ind(size(cell_edge_adj), edgelist(:, 2)', 1:N_edge)) = 1;

    % Compute vertex IDs defining each edge
    edge_vid = zeros(N_edge, 2); % Preallocate edge-vertex matrix
    for i_e = 1:N_edge % Loop through all edges
        edge_vid(i_e, :) = find(all(cell_vertex_adj(cell_edge_adj(:, i_e) > 0, :), 1));
    end
end

% If edge vectors are needed, compute their lengths and directions
if nargout > 3 && nargin > 1
    edge_vector = zeros(N_edge, 2); % Preallocate vector list
    for i_e = 1:N_edge % Loop through edges
        two_vid = edge_vid(i_e, :); % Get vertices defining the edge
        two_end_position = pbc_relocate(vertex_position(two_vid(1), :), ...
                                        vertex_position(two_vid, :), box_size); % Handle periodic boundaries
        edge_vector(i_e, :) = two_end_position(1, :) - two_end_position(2, :); % Compute vector
    end
end

% If edge tensions are needed, compute them
if nargout > 4 && nargin > 3
    % Compute the perimeter of each cell using edge lengths
    cell_perimeter_list = cell_edge_adj * vecnorm(edge_vector, 2, 2);
    % Compute deviation of cell perimeters from preferred values
    cell_peri_devia_list = 2 * K_P .* (cell_perimeter_list - P0_list);
    % Calculate tension for each edge
    edge_tension_list = cell_edge_adj.' * cell_peri_devia_list;
end

end % End of function

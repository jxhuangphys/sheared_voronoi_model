function [num_t1, losing_pair, gaining_pair] = ...
    get_happen_t1(old_adj, new_adj, cell_ind)
% Function to detect T1 transitions (neighbor exchanges) in a tissue
% Inputs:
%   old_adj - Adjacency matrix before the transition
%   new_adj - Adjacency matrix after the transition
%   cell_ind (optional) - Indices of cells to consider for transitions
% Outputs:
%   num_t1 - Number of T1 transitions (losing pairs)
%   losing_pair - List of cell pairs that lost adjacency
%   gaining_pair (optional) - List of cell pairs that gained adjacency

% If cell indices are provided, reduce the adjacency matrices to those cells
if nargin == 3
    cell_ind = sort(cell_ind(:)).'; % Ensure cell_ind is a row vector
    old_adj = old_adj(cell_ind, cell_ind); % Subset the old adjacency matrix
    new_adj = new_adj(cell_ind, cell_ind); % Subset the new adjacency matrix
end

% Size of the adjacency matrix
mat_size = size(old_adj, 1);

% --- Detect losing pairs ---
% Find indices where adjacency exists in old_adj but not in new_adj
losing_ids = find(old_adj - new_adj == 1); 
% Convert linear indices to row and column subscripts
[losing_1, losing_2] = ind2sub([mat_size mat_size], losing_ids);

% Remove duplicates by ensuring pairs are sorted as [small, big]
unique_ind = find(losing_1 - losing_2 < 0);
losing_1 = losing_1(unique_ind);
losing_2 = losing_2(unique_ind);

% Combine losing pairs into a single array
losing_pair = [losing_1 losing_2];
% Count the number of losing pairs
num_t1 = size(losing_pair, 1);

% If specific cell indices were used, map back to the global indices
if exist('cell_ind', 'var')
    losing_pair = cell_ind(losing_pair);
end

% --- Detect gaining pairs (if requested) ---
if nargout == 3
    % Find indices where adjacency does not exist in old_adj but exists in new_adj
    gaining_ids = find(old_adj - new_adj == -1);
    % Convert linear indices to row and column subscripts
    [gaining_1, gaining_2] = ind2sub([mat_size mat_size], gaining_ids);

    % Remove duplicates by ensuring pairs are sorted as [small, big]
    unique_ind = find(gaining_1 - gaining_2 < 0);
    gaining_1 = gaining_1(unique_ind);
    gaining_2 = gaining_2(unique_ind);

    % Combine gaining pairs into a single array
    gaining_pair = [gaining_1 gaining_2];

    % If specific cell indices were used, map back to the global indices
    if exist('cell_ind', 'var')
        gaining_pair = cell_ind(gaining_pair);
    end
end

end % End of the function

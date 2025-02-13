function draw_shear_voronoi(center_xy, cell_chain, vertex_position, gamma, box_size)
    % Function to draw shear-transformed Voronoi diagrams with periodic boundary conditions
    %
    % <input>
    % center_xy: N_cell x 2 array containing the coordinates of the Voronoi cell centers
    % cell_chain: Cell array, where each element contains indices of vertices forming a cell
    % vertex_position: N_vertex x 2 array containing the coordinates of vertices
    % gamma: Shear factor for the Lees-Edwards boundary condition
    % box_size: [Lx, Ly], the dimensions of the periodic boundary box
    
    % Get the number of Voronoi cells
    N_cell = numel(cell_chain);

    % Set up the axis properties for visualization
    axis equal  % Ensure equal scaling for x and y axes
    hold on     % Keep all plots on the same figure

    % Draw the periodic boundary box
    rectangle('position', [0 0 box_size(1) box_size(end)]) % Box boundaries
    plot(center_xy(:, 1), center_xy(:, 2), 'r.')          % Plot cell centers as red dots

    % Iterate through each cell to draw its Voronoi polygon
    for i_c = 1:N_cell
        % Get the vertex positions for the current cell
        chain_xy = vertex_position(cell_chain{i_c}, :);

        % Apply the Lees-Edwards boundary condition for shear
        % Adjust vertices that cross the shear boundaries based on their y-coordinates
        if center_xy(i_c, 2) > box_size * 2/3
            % Cells in the upper boundary: apply positive shear adjustment
            chain_xy(chain_xy(:, 2) <= box_size / 3, 1) = ...
                chain_xy(chain_xy(:, 2) <= box_size / 3, 1) + gamma * box_size;
        elseif center_xy(i_c, 2) < box_size / 3
            % Cells in the lower boundary: apply negative shear adjustment
            chain_xy(chain_xy(:, 2) >= box_size * 2/3, 1) = ...
                chain_xy(chain_xy(:, 2) >= box_size * 2/3, 1) - gamma * box_size;
        end

        % Relocate vertices to the primary periodic box using a helper function
        chain_xy = pbc_relocate(center_xy(i_c, :), chain_xy, box_size);

        % Plot the closed polygon representing the Voronoi cell
        plot([chain_xy(:, 1); chain_xy(1, 1)], [chain_xy(:, 2); chain_xy(1, 2)], ...
            'color', 'b', 'linewidth', 1); % Blue lines for cell borders
    end
end

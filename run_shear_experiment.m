% Parameters defining the tissue and its energy model
p0 = 3.78; % Target shape index for cells
N_cell = 64; % Number of cells in the simulation
box_size = sqrt(N_cell); % Size of the simulation box

ka = 0; % Area elasticity constant for cells

% Relaxation parameters for the minimization process
force_tol = 1e-13; % Tolerance for force convergence during relaxation
fire_step = 1.5e4; % Number of steps for the FIRE relaxation algorithm
K_A = ka * ones(N_cell, 1); % Uniform area elasticity for all cells
K_P = 1 * ones(N_cell, 1); % Uniform perimeter elasticity for all cells
A0_list = 1 * ones(N_cell, 1); % Preferred area for all cells
P0_list = p0 * ones(N_cell, 1); % Preferred perimeter for all cells

% Shear parameters
final_gamma = 6; % Total shear strain to be applied
shear_step = 3000; % Number of steps for incremental shear application
d_gamma = final_gamma / shear_step; % Increment of shear strain per step

% Random seed for reproducibility
seed_idx = 1; % Seed index
seed = RandStream('mcg16807', 'Seed', seed_idx); % Initialize random stream
RandStream.setGlobalStream(seed); % Set global stream for random number generation

% Initialize random positions for cell centers within the box
center_xy = rand(N_cell, 2) * box_size;

% Relax to a ground state using FIRE algorithm
[center_xy, cell_chain, vertex_position, average_energy, ...
 mean_force, area_list, perimeter_list] = ...
 minimize_shear_voronoi_FIRE(center_xy, 0, K_A, A0_list, K_P, P0_list, ...
 3 * fire_step, force_tol);

gamma = 0; % Initial shear strain

data = struct('center_xy',[],'cell_chain',[],'vertex_position',[],'current_strain',[]);% structure that store the simulation data

% Main loop to apply incremental shear and record results
for i_s = 1:shear_step + 1
    % Relax cell positions using FIRE after applying shear
    [center_xy, cell_chain, vertex_position, average_energy, ...
     mean_force, area_list, perimeter_list] = ...
     minimize_shear_voronoi_FIRE(center_xy, gamma, K_A, A0_list, K_P, P0_list, ...
     fire_step, force_tol);

     % store data
    data(i_s).center_xy = center_xy;
    data(i_s).cell_chain = cell_chain;
    data(i_s).vertex_position = vertex_position;
    data(i_s).gamma = gamma;

	% Draw the current state of the Voronoi tessellation
    figure; % Create a new figure for each state
    draw_shear_voronoi(center_xy, cell_chain, vertex_position, gamma, box_size);
    
    % Save the figure as a PNG image
    saveas(gcf, ['qs_ka', num2str(ka, 4), '_N', num2str(N_cell), '_p', num2str(p0), ...
                 '_gamma', num2str(gamma, 4), '.png']);
    close; % Close the figure to save memory

    
    % Update shear strain and apply simple shear to cell positions
    gamma = gamma + d_gamma;
    center_xy(:, 1) = center_xy(:, 1) + d_gamma * (center_xy(:, 2) - box_size / 2); % Apply shear
    center_xy = mod(center_xy, box_size); % Periodic boundary conditions

    % Update Voronoi tessellation to account for shear
    [cell_chain, vertex_position, area_list, perimeter_list] = ...
        make_point_voronoi_shear_lebc(center_xy, gamma);
end

%%%%% saving the data
save( ['qs_ka', num2str(ka, 4), '_N', num2str(N_cell), '_p', num2str(p0),'.mat'],'data')




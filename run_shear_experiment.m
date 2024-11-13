% Parameters defining the tissue and its energy model
p0 = 0.78; % Target shape index for cells
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

% Structures to store simulation parameters and results
qs_parameter = struct(); % Parameters for quasi-static shear simulations
qs_shear_result = struct(); % Results of the shear simulations

% Shear parameters
final_gamma = 6; % Total shear strain to be applied
shear_step = 3000; % Number of steps for incremental shear application
d_gamma = final_gamma / shear_step; % Increment of shear strain per step

% Assigning parameters to the parameter structure
qs_parameter.N_cell = N_cell;
qs_parameter.p0 = p0;
qs_parameter.d_gamma = d_gamma;
qs_parameter.force_tol = force_tol;
qs_parameter.shear_step = shear_step;
qs_parameter.final_gamma = final_gamma;
qs_parameter.fire_step = fire_step;

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

% Record adjacency information of the cells
old_adj = get_tissue_adj(cell_chain);

% Data structures to record shear results over steps
snapshot = struct(); % Snapshot of simulation state
gamma_list = zeros(shear_step + 1, 1); % List of applied shear strains
stress_int = zeros(shear_step + 1, 3); % Shear stress components
shape_tensor_list = zeros(shear_step + 1, 3); % Global shape tensor components
shape_ind_mean_std = zeros(shear_step + 1, 2); % Mean and standard deviation of shape indices
t1_history = zeros(shear_step + 1, 1); % T1 transitions history
long_axis_info = zeros(shear_step + 1, 6); % Information about cell shapes and orientation
mean_force_list = zeros(shear_step + 1, 1); % Mean forces at each step
gamma = 0; % Initial shear strain
total_t1 = 0; % Total number of T1 transitions
start_step = 1; % Start step for simulation

% Main loop to apply incremental shear and record results
for i_s = start_step:shear_step + 1
    % Relax cell positions using FIRE after applying shear
    [center_xy, cell_chain, vertex_position, average_energy, ...
     mean_force, area_list, perimeter_list] = ...
     minimize_shear_voronoi_FIRE(center_xy, gamma, K_A, A0_list, K_P, P0_list, ...
     fire_step, force_tol);

    % Detect T1 transitions (neighbor exchanges)
    new_adj = get_tissue_adj(cell_chain);
    [happen_t1, losing_pair] = get_happen_t1(old_adj, new_adj);
    total_t1 = total_t1 + happen_t1;
    old_adj = new_adj;

    % Calculate stress and shape parameters
    mean_stress_int = get_shear_tissue_stress(center_xy, cell_chain, vertex_position, ...
                    K_A, A0_list, K_P, P0_list, gamma, box_size);

    [centroid_xy, long_axis_angle, global_long_angle, ...
     order_parameter, aspect_ratio, global_shape_tensor] = ...
     get_shear_voronoi_axis(center_xy, cell_chain, ...
     vertex_position, gamma, box_size);

    % Record simulation state
    snapshot(i_s).center_xy = center_xy;
    snapshot(i_s).losing_pair = losing_pair;

    % Store calculated data for analysis
    gamma_list(i_s, 1) = gamma;
    mean_force_list(i_s, 1) = mean_force;
    stress_int(i_s, :) = mean_stress_int([1 2 4]);
    shape_tensor_list(i_s, :) = global_shape_tensor([1 2 4]);
    shape_index = perimeter_list ./ sqrt(area_list);
    shape_ind_mean_std(i_s, :) = [mean(shape_index), std(shape_index)];
    t1_history(i_s) = total_t1;
    shape_anisotropy = (aspect_ratio - 1) ./ (aspect_ratio + 1);
    long_axis_info(i_s, :) = [global_long_angle, order_parameter, ...
        mean(aspect_ratio), median(aspect_ratio), ...
        mean(shape_anisotropy), median(shape_anisotropy)];

    % Update shear strain and apply simple shear to cell positions
    gamma = gamma + d_gamma;
    center_xy(:, 1) = center_xy(:, 1) + d_gamma * (center_xy(:, 2) - box_size / 2); % Apply shear
    center_xy = mod(center_xy, box_size); % Periodic boundary conditions

    % Update Voronoi tessellation to account for shear
    [cell_chain, vertex_position, area_list, perimeter_list] = ...
        make_point_voronoi_shear_lebc(center_xy, gamma);
end

% Save results and parameters to file
qs_shear_result(seed_idx).seed = rng;
qs_shear_result(seed_idx).current_gamma = gamma;
qs_shear_result(seed_idx).gamma_list = gamma_list;
qs_shear_result(seed_idx).mean_force_list = mean_force_list;
qs_shear_result(seed_idx).t1_history = t1_history;
qs_shear_result(seed_idx).long_axis = long_axis_info;
qs_shear_result(seed_idx).stress_int = stress_int;
qs_shear_result(seed_idx).shape_tensor_list = shape_tensor_list;
qs_shear_result(seed_idx).shape_index = shape_ind_mean_std;
qs_shear_result(seed_idx).snapshot = snapshot;

file_data(['fire/qs_ka', num2str(ka, 4), '_N', num2str(N_cell), '_p', num2str(p0), ...
           '_gamma', num2str(gamma), ...
           '_dgamma', num2str(d_gamma), ...
           '_relaxstep', num2str(fire_step), '.mat'], ...
           qs_shear_result, qs_parameter);

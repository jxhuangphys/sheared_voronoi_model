function [center_xy, cell_chain, vertex_position, ...
          average_energy, mean_force, area_list, perimeter_list] = ...
          minimize_shear_voronoi_FIRE(center_xy, gamma, K_A, ...
          A0_list, K_P, P0_list, max_step, force_tol)

% Initialize the number of cells and the simulation box size
N_cell = size(center_xy, 1); % Number of cells
box_size = sqrt(N_cell); % Size of the simulation box

% Ensure gamma (shear strain) is in the range [-0.5, 0.5] for periodic boundaries
gamma = mod(gamma, 1);
if gamma > 0.5 
    gamma = gamma - 1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FIRE Algorithm Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
threshold_step = 5; % Steps before increasing the time step
finc = 1.1; % Multiplicative factor to increase time step
fdec = 0.5; % Multiplicative factor to decrease time step
f_alpha = 0.99; % Decay factor for the FIRE damping parameter
alpha_start = 0.1; % Initial FIRE damping parameter
dt_max = 0.1; % Maximum time step
dt_min = 0.001; % Minimum time step
fire_mass = 4; % Artificial mass for FIRE velocity updates

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialization of FIRE algorithm variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dt = sqrt(dt_max * dt_min); % Initial time step
center_velocity = zeros(N_cell, 2); % Initialize velocities of cell centers
alpha = alpha_start; % Start with initial FIRE damping
mean_force = 1e8; % Initialize mean force to a high value
loop_len = 100; % Length of each inner loop
num_loop = ceil(max_step / loop_len); % Number of outer loops
right_step = 0; % Counter for successful steps

% Main outer loop for FIRE relaxation
for i_loop = 1:num_loop
    
    for i_s = 1:loop_len
        % Compute Voronoi tessellation and cell properties under shear
        [cell_chain, vertex_position, area_list, perimeter_list] = ...
                        make_point_voronoi_shear_lebc(center_xy, gamma);

        % Compute forces on cell centers based on energy gradients
        center_force = get_shear_voronoi_force(center_xy, cell_chain, ...
                        vertex_position, gamma, ...
                        K_A, A0_list, K_P, P0_list, area_list, perimeter_list);

        % Normalize forces to compute the direction
        force_hat = center_force ./ vecnorm(center_force, 2, 2);
        force_hat(isnan(force_hat)) = 0; % Handle NaN cases from zero force

        % Compute the power (dot product of force and velocity)
        Pfire = sum(sum(center_velocity .* center_force));

        % Adjust time step and velocities based on FIRE conditions
        if (Pfire > 0 && right_step > threshold_step) || ...
           (i_loop > num_loop/2 && dt < mean_force / 100)
            % Increase time step and reduce damping
            dt = min(dt * finc, dt_max);
            alpha = alpha * f_alpha;
            vel_mag = vecnorm(center_velocity, 2, 2);
            center_velocity = (1 - alpha) * center_velocity + alpha * vel_mag .* force_hat;

        elseif Pfire <= 0
            % Reset velocities and reduce time step
            dt = max(dt_min, dt * fdec);
            alpha = alpha_start;
            center_velocity = zeros(N_cell, 2); % Reset velocity
            right_step = 0; % Reset step counter
        end
        
        % Update velocity and position of cells
        right_step = right_step + 1;
        center_velocity = center_velocity + center_force / fire_mass * dt;
        center_xy = center_xy + center_velocity * dt;

        % Apply periodic boundary conditions with shear deformation
        center_xy(:, 1) = center_xy(:, 1) - box_size * gamma * ...
                                    ((center_xy(:, 2) > box_size) - (center_xy(:, 2) < 0));
        center_xy = mod(center_xy, box_size); % Wrap positions into the box

        % Compute mean force to check for convergence
        mean_force = mean(vecnorm(center_force, 2, 2));
    end
    
    % Compute average energy of the system
    average_energy = mean(K_P .* (perimeter_list - P0_list).^2 ...
            + K_A .* (area_list - A0_list).^2);

    % Check minimum cell distance for potential instabilities
    min_dist = min(pdist(center_xy));
        
    % Terminate if mean force is below the threshold
    if mean_force < force_tol
        break
    end
    
    % Adjust time step if cells overlap excessively
    if dt == dt_min && min_dist < 0.25
        dt = dt * 10;
    end
    
    % Gradually increase the maximum time step to speed up convergence
    if right_step > max_step / 10
        dt_max = dt_max * 1.05;
        right_step = threshold_step + 1;
    end

end

end

function points = pbc_relocate(center, points, box_size)
% Function to relocate points to their periodic boundary conditions (PBC)
% relative to a given center. Ensures that points are confined within a
% simulation box defined by box_size and are close to the center.

% Iterate over each point in the points array
for i = 1:size(points, 1)
    % Extract the current point coordinates
    i_point = points(i, :);
    
    % Calculate the relative difference between the point and the center
    % Normalize this difference by dividing by the box size and shift by 0.5
    % to center it in the range [0, 1].
    difference = (i_point - center)./box_size + 0.5;
    
    % Adjust the point coordinates to fall within the box boundaries.
    % By taking the floor of the normalized difference, we determine the
    % number of box lengths to subtract from the point.
    % This ensures the final point remains within the range [0, box_size].
    points(i, :) = i_point - floor(difference).*box_size;
end

% End of function
end

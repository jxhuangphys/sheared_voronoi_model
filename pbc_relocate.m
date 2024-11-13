function points = pbc_relocate(center, points, box_size)

for i = 1:size(points, 1)
    i_point = points(i, :);
    difference = (i_point - center)./box_size + 0.5;
     
    % we want this difference to be between 0 and 1
    % because (difference - 0.5) should be between -0.5 and 0.5
    points(i, :) = i_point - floor(difference).*box_size;
end
 

end % end of function
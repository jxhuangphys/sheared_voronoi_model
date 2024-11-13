function [centroid_xy, long_axis_angle, global_long_angle, ...
    nematic_order, aspect_ratio, shape_tensor] = ...
                            get_shear_voronoi_axis(center_xy, cell_chain, vertex_position, ...
                            gamma, box_size)
                       
                        
N_cell = numel(cell_chain);
if ~exist('box_size', 'var')
    box_size = sqrt(N_cell);
end
% Lx = box_size(1);
Ly = box_size(end);

shape_tensor = zeros(N_cell, 3);
centroid_xy = zeros(N_cell, 2);
long_axis_angle = zeros(N_cell, 1);
aspect_ratio = zeros(N_cell, 1);
for i_c = 1:N_cell
    chain_xy = vertex_position(cell_chain{i_c}, :);
    if gamma && center_xy(i_c, 2)>Ly*2/3
        chain_xy(chain_xy(:, 2)<=Ly/3, 1) = chain_xy(chain_xy(:, 2)<=Ly/3, 1) + gamma*Ly;
    elseif gamma && center_xy(i_c, 2)<Ly/3
        chain_xy(chain_xy(:, 2)>=Ly*2/3, 1) = chain_xy(chain_xy(:, 2)>=Ly*2/3, 1) - gamma*Ly;
    end
    chain_xy = pbc_relocate(center_xy(i_c, :), chain_xy, box_size);

    [ geom, iner, cpmo ] = polygeom( chain_xy(:, 1), chain_xy(:, 2) );
    
    centroid_xy(i_c, :) = geom(2:3);
    long_axis_angle(i_c) = mod(cpmo(2), pi);
    
    shape_tensor(i_c, 1) = iner(5); 
    shape_tensor(i_c, 2) = iner(6);
    shape_tensor(i_c, 3) = iner(4);

    aspect_ratio(i_c) = sqrt(cpmo(3)/cpmo(1));
end


global_angle_tensor = [mean(cos(long_axis_angle).^2) mean(cos(long_axis_angle).*sin(long_axis_angle)); ...
    mean(cos(long_axis_angle).*sin(long_axis_angle)) mean(sin(long_axis_angle).^2)];
[eigen_vec, ~] = eig(global_angle_tensor);
global_long_axis = eigen_vec(:, 2);
global_long_angle = atan(global_long_axis(2)/global_long_axis(1));
nematic_order = mean(cos(2*(long_axis_angle-global_long_angle)));


    function [ geom, iner, cpmo ] = polygeom( x, y ) 
    %   POLYGEOM Geometry of a planar polygon
    %   reference : 
    %   https://paperpool.github.io/On_the_Calculation_of_Moments_of_Polygons.pdf
    %
    %   POLYGEOM( X, Y ) returns area, X centroid,
    %   Y centroid and perimeter for the planar polygon
    %   specified by vertices in vectors X and Y.
    %
    %   [ GEOM, INER, CPMO ] = POLYGEOM( X, Y ) returns
    %   area, centroid, perimeter and area moments of 
    %   inertia for the polygon.
    %   GEOM = [ area   X_cen  Y_cen  perimeter ]
    %   INER = [ Ixx    Iyy    Ixy    Iuu    Ivv    Iuv ]
    %     u,v are centroidal axes parallel to x,y axes.
    %   CPMO = [ I1     ang1   I2     ang2   J ]
    %     I1,I2 are centroidal principal moments about axes
    %         at angles ang1,ang2.
    %     ang1 and ang2 are in radians.
    %     J is centroidal polar moment.  J = I1 + I2 = Iuu + Ivv

    
    
    % begin function POLYGEOM

    % check if inputs are the same size
    if ~isequal( size(x), size(y) )
      error( 'X and Y must be the same size');
    end

    % temporarily shift data to mean of vertices for improved accuracy
    xm = mean(x);
    ym = mean(y);
    x = x - xm;
    y = y - ym;

    % summations for counter-clockwise boundary
    xp = x( [2:end 1] );
    yp = y( [2:end 1] );
    a = x.*yp - xp.*y;

    A = sum( a ) /2;
    xc = sum( (x+xp).*a  ) /6/A;
    yc = sum( (y+yp).*a  ) /6/A;
    Ixx = sum( (y.*y +y.*yp + yp.*yp).*a  ) /12;
    Iyy = sum( (x.*x +x.*xp + xp.*xp).*a  ) /12;
    Ixy = sum( (x.*yp +2*x.*y +2*xp.*yp + xp.*y).*a  ) /24;

    dx = xp - x;
    dy = yp - y;
    P = sum( sqrt( dx.*dx +dy.*dy ) );

    % check for CCW versus CW boundary
    if A < 0
      A = -A;
      Ixx = -Ixx;
      Iyy = -Iyy;
      Ixy = -Ixy;
    end

    % centroidal moments
    Iuu = Ixx - A*yc*yc;
    Ivv = Iyy - A*xc*xc;
    Iuv = Ixy - A*xc*yc;
    J = Iuu + Ivv;

    % replace mean of vertices
    x_cen = xc + xm;
    y_cen = yc + ym;
    Ixx = Iuu + A*y_cen*y_cen;
    Iyy = Ivv + A*x_cen*x_cen;
    Ixy = Iuv + A*x_cen*y_cen;

    % principal moments and orientation
    I = [ Iuu  -Iuv ;
         -Iuv   Ivv ];
    [ eig_vec, eig_val ] = eig(I);
    I1 = eig_val(1, 1);
    I2 = eig_val(2, 2); % I2 > I1
    ang1 = atan2( eig_vec(2,1), eig_vec(1,1) );
    ang2 = atan2( eig_vec(2,2), eig_vec(1,2) );

    % return values
    geom = [ A  x_cen  y_cen  P ];
    iner = [ Ixx  Iyy  Ixy  Iuu  Ivv  Iuv ];
    cpmo = [ I1  ang1  I2  ang2  J ];

    % bottom of polygeom

    end % end of polygeom



end % end of the whole function



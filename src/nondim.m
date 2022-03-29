function [x_nondim, A, DU, VU] = nondim(x_dim, frame)
%NONDIM non-dimensionalize a state vector in the ECI or SCI frame
% 
% [x_nondim, J, DU, VU] = NONDIM(x_dim, frame)
% 
% Inputs:   x_dim [km; km/s] (6xN) dimensional state vector
%           frame [] (string) "ECI" or "SCI"
% 
% Outputs:  x_nondim [] (6x1) non-dimensional state vector
%           A [] (6x6) transformation matrix from dimensional to non-dimensional
%           DU [km] (scalar) distance unit
%           VU [km/s] (scalar) velocity unit
% 
% See also: 

% Author: Jared Blanchard 	Date: 2022/02/23 13:20:30 	Revision: 0.1 $

x_dim = make_column_vecs(x_dim);

if frame == "ECI"
    GM = cspice_bodvrd( 'EARTH', 'GM', 1 );
    R = cspice_bodvrd( 'EARTH', 'RADII', 3); %[km] radii of Earth Ellipsoid
    DU = R(1); % distance unit is radius of the earth
elseif frame == "SCI"
    GM = cspice_bodvrd( 'SUN', 'GM', 1 );
    AU = cspice_convrt(1, 'AU', 'km'); %[km]
    DU = AU; % distance unit is astronomical unit
else
    error("frame should be 'ECI' or 'SCI'")
end

VU = sqrt(GM/DU); % velocity unit
A = [eye(3)/DU,  zeros(3);
      zeros(3), eye(3)/VU];
x_nondim = A*x_dim;

end

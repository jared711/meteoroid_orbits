function [r] = azel2cart(az, el, d)
%AZEL2R converts azimuth, elevation, and magnitude measurements into a
%cartesian position vector with the same units as the magnitude
% 
% [r] = AZEL2CART(az, el, d)
% 
% Inputs:   az = [deg] (scalar) azimuth
%           el = [deg] (scalar) elevation angle
%           d  = [km]  (scalar) magnitude {1} 
% 
% Outputs:  r  = [km]  (3x1)    position vector
% 
% See also:

% Author: Jared Blanchard 	Date: 2021/01/15 13:54:15 	Revision: 0.1 $
% Author: Jared Blanchard 	Date: 2021/03/2 14:34:23 	Revision: 0.2 $

if nargin < 3;  d = 1;              end % if no magnitude specified, make a unit vector


r = [d*cosd(el)*sind(az);
     d*cosd(el)*cosd(az);
     d*sind(el)];

end

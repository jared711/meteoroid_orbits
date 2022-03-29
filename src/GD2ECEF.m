function r_ECEF = GD2ECEF(lat,lon,alt)
%GD2ECEF converts Geodetic latitude (deg), longitude (deg),
% and altitude to cartesian position vector in ECEF frame
% 
% [r_ECEF] = GD2ECEF(lat, lon, alt)
% 
% Inputs:   lat = [rad/deg] (scalar) geodetic latitude
%           lon = [rad/deg] (scalar) longitude
%           alt = [km]  (scalar) altitude {0}
% 
% Outputs:  r_ECEF = [km] (3x1) position vector in ECEF frame
% 
% See also: GC2ECEF 

% Author: Jared Blanchard 	Date: 2021/01/18 21:01:32 	Revision: 0.1 $

if nargin < 3;  alt = 0;            end

% these values are take from the WGS84 ellipsoid
% wgs84Ellipsoid('km').SemimajorAxis
% wgs84Ellipsoid('km').Eccentricity
r_E = 6378.137; % [km] equatorial radius of the Earth
e_E = 0.081819190842621;   % [] eccentricity of the Earth ellipsoid

N = r_E/sqrt(1-(e_E^2)*sind(lat)^2);

r_ECEF = [(N + alt)*cosd(lat)*cosd(lon);
          (N + alt)*cosd(lat)*sind(lon);
          (N*(1-e_E^2) + alt)*sind(lat)];

end

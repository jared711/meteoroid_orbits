function [p] = p_srp(d)
%P_SRP solar radiation pressure at distance d from the sun
% 
% [p] = P_SRP(d)
% 
% Inputs:   d [km] distance from sun
% 
% Outputs:  p [N/m^2] solar radiation pressure
% 
% See also: 

% Author: Jared Blanchard 	Date: 2021/03/12 15:28:50 	Revision: 0.1 $


G_sc = 1360.8; % [W/m^2] Solar light flux measured at 1 AU
c = 2.99792458e8; % [m/s] speed of light
AU = cspice_convrt(1, 'AU', 'km'); %[km]
p = (G_sc*(AU)^2)/(c*d^2); % [Pa or kg/m/s^2]

end

% Constants from Kopp, G.; Lean, J. L. (2011). "A new, lower value of total solar irradiance: Evidence and climate significance". Geophysical Research Letters. 38 (1): n/a. Bibcode:2011GeoRL..38.1706K. doi:10.1029/2010GL045777.

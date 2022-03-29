function [x_SCI, R] = ECI2SCI(x_ECI, et)
%ECI2SCI converts a cartesian state vector from ECI to SCI frames and
%outputs the corresponding rotation matrix
% 
% [x_SCI, R] = ECI2SCI(x_ECI, et)
% 
% Inputs:   x_ECI = [km;km/s] (6x1) full state vector of meteor in ECI frame 
%           et = (scalar) ephemeris time
% 
% Outputs:  x_SCI = [km;km/s] (6x1) full state in SCI frame
% 
% See also: 

% Author: Jared Blanchard 	Date: 2022/01/28 12:36:57 	Revision: 0.1 $


% Rotation matrix from ECI to SCI
if ~iscolumn(x_ECI);    x_ECI = x_ECI'; end

R = cspice_sxform('J2000','ECLIPJ2000',et); % transform from J2000 to ECLIPJ2000 at epoch et

x_sun_ECI = cspice_spkezr('SUN', et, 'J2000', 'NONE', 'EARTH'); % [km; km/s] state of sun in ECI frame

x_SCI = R*(x_ECI - x_sun_ECI);

end

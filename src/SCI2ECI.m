function [x_ECI, R] = SCI2ECI(x_SCI, et)
%ECI2SCI converts a cartesian state vector from SCI to ECI frames and
%outputs the corresponding rotation matrix
% 
% [x_ECI, R] = ECI2SCI(x_SCI, et)
% 
% Inputs:   x_SCI = [km;km/s] (6x1) full state vector of meteor in ECI frame 
%           et = (scalar) ephemeris time
% 
% Outputs:  x_ECI = [km;km/s] (6x1) full state in SCI frame
% 
% See also: 

% Author: Jared Blanchard 	Date: 2022/01/28 12:36:57 	Revision: 0.1 $

% Rotation matrix from SCI to ECI
R = cspice_sxform('ECLIPJ2000','J2000',et); % transform from ECLIPJ2000 to J2000 at epoch et

x_sun_ECI = cspice_spkezr('SUN', et, 'J2000', 'NONE', 'EARTH'); % [km; km/s] state of sun in ECI frame

x_ECI = R*x_SCI + x_sun_ECI;
end
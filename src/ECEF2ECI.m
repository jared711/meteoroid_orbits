 function [x_ECI, R] = ECEF2ECI(x_ECEF, time)
%ECEF2ECI converts a cartesian state vector from ECEF frame to ECI frame
%given a time vector [YYYY;MM;DD;hh;mm;ss]
% 
% [x_ECI, R] = ECEF2ECI(x_ECEF, time)
% 
% Inputs:   x_ECEF = [km; km/s] (6x1) state in ECEF frame 
%           time = [YYYY;MM;DD;hh;mm;ss] (6x1) time vector
% 
% Outputs:  rv_ECI = [km; km/s] (6x1) state in ECI frame
%           R = [] (6x6) rotation matrix from ECEF to ECI frame
% 
% Author: Jared Blanchard 	Date: 2021/01/15 15:39:59 	Revision: 0.1 $

if isrow(x_ECEF);       x_ECEF = x_ECEF';                       end
if ~iscolumn(x_ECEF);   error('rv_ECEF should be a vector');    end

% ephemeris time conversion from SPICE
et = cspice_str2et(datestr(time)); % [sec] ephemeris time

% rotation matrix from ECEF to ECI frame
R = cspice_sxform('IAU_EARTH','J2000',et);

x_ECI = R*x_ECEF;

end
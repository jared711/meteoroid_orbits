function [x_SCI, S] = ENU2SCI(x_ENU, phi, lambda, h, time)
%ENU2SCI Converts a 6D cartesian state vector from the ENU frame to the SCI
%frame and outputs the jacobian matrix
% 
% [x_SCI, S] = ENU2SCI(x_ENU, phi, lambda, h, time)
% 
% Inputs:   x_ENU {6x1} [km;km/s] cartesian state vector in ENU frame
%           phi = [deg] (scalar) latitude of observer radar
%           lambda = [deg] (scalar) longitude of observer radar
%           h = [km] (scalar) altitude of observer radar above geoid
%           time = [YYYY;MM;DD;hh;mm;ss] (6x1) time vector
% 
% Outputs:  x_SCI = [km;km/s] (6x1) cartesian state vector in SCI frame
% 
% See also: 

% Author: Jared Blanchard 	Date: 2022/01/28 11:37:47 	Revision: 0.1 $

% ephemeris time conversion from SPICE
et = cspice_str2et(datestr(time)); % [sec] ephemeris time

[x_ECEF, R1] = ENU2ECEF(x_ENU, phi, lambda, h);
[x_ECI, R2] = ECEF2ECI(x_ECEF, et);
[x_SCI, R3] = ECI2SCI(x_ECI, et);
S = R3*R2*R1;

end

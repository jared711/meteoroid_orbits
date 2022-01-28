function [x_ECEF, R] = ENU2ECEF(x_ENU, phi, lambda, h)
%ENU2ECEF converts a cartesian state vector from the ENU frame located at 
% lat,lon in rad/deg) into the ECEF frame
% 
% [x_E] = ENU2ECEF(x_ENU, lat, lon, h, geodetic, ang_unit)
% 
% Inputs:   x_ENU = [km;km/s] (6x1) full state in ENU frame
%           lat = [deg] (scalar) latitude of observer
%           lon = [deg] (scalar) longitude of observer
%           h = [km] (scalar) height of observer above geoid
% 
% Outputs:  x_ECEF = [km;km/s] (6x1) full state in ECEF frame
%           R = [] (6x6) rotation matrix from ENU to ECEF frame
% 
% See also: 

% Author: Jared Blanchard 	Date: 2022/01/28 11:50:05 	Revision: 0.1 $

if nargin < 4;      h = 0;              end
if isrow(x_ENU);    x_ENU = x_ENU';     end
if ~iscolumn(x_ENU);    error('x_ENU should be a vector');  end

% A = radar_rot(lat, lon); % [] rotation matrix from R to E
RzRx = rotz(pi/2 + lambda)*rotx(pi/2 - phi);
R = [    RzRx, zeros(3);
     zeros(3),     RzRx];

% geodetic radar location in ECEF cartesian using SPICE
radii = cspice_bodvrd( 'EARTH', 'RADII', 3 );
flat = (radii(1) - radii(3))/radii(1);
r_radar_ECEF = cspice_georec(lambda, phi, h, radii(1), flat);

x_radar_ECEF = [r_radar_ECEF; zeros(3,1)];

x_ECEF = R*x_ENU + x_radar_ECEF;

end

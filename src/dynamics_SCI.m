function [vec_dot] = dynamics_SCI(t, vec, et0, flags, A, m, bodies)
%DYNAMICS_SCI calculates the time derivative of a
%state vector given in cartesian coordinates in an inertial frame about the
%central body. Can take drag and solar radiation pressure into account.
% 
% [vec_dot] = dynamics_SCI(t, vec, flags, A, m)
% 
% Inputs:   t [sec] (scalar) time used by ode45, ode113 etc...
%           vec [km;km/s;NON] (42x1) cartesian state vector
%           flags [] (3x1) flags = [flag_drag; flag_SRP; flag_TB];
%           A = [m^2] (scalar) cross-sectional area
%           m = [kg] (scalar) mass of meteoroid
%           options (struct)
%               central_body (struct) central gravitational body
%               meteor (struct) object/spacecraft
%               flag_j2 (bool) indicator of j2 {0}
%               flag_3BP (bool) indicator of third body {0}
%               flag_drag (bool) indicator that we are using atmospheric drag {0}
%               flag_srp (bool) indicator that we are using solar radiation pressure {0}
%               flag_poynt (bool) indicator of poynting-robertson effect {0}
%               flag_YORP (bool) indicator of Yarkovsky ORP effect {0}
%               flag_lorentz (bool) indicator of lorentz forces {0}
% 
% Outputs:  x_dot [km/s;km/s^2] (6x1) d/dt of cartesian state vector
% 
% See also: 

% Author: Jared Blanchard 	Date: 2022/01/28 13:58:56 	Revision: 0.1 $

if isrow(vec);    vec = vec';                     end

et = et0 + t; % [sec] ephemeris time
gm = cspice_bodvrd( 'SUN', 'GM', 1 );

x_SCI = vec(1:6); % [km;km/s] state vector
r_SCI = vec(1:3); % [km] position vector
v_SCI = vec(4:6); % [km/s] velocity vector
a_SCI = (-gm/(norm(r_SCI)^3))*r_SCI; % [km/s^2] acceleration vector

PHI = reshape(vec(7:42),6,6);

d_drdot = zeros(3);
d_dr = gm/(norm(r_SCI)^5)*...
   [3*r_SCI(1)^2 - norm(r_SCI)^2, 3*r_SCI(1)*r_SCI(2), 3*r_SCI(1)*r_SCI(3);
    3*r_SCI(1)*r_SCI(2), 3*r_SCI(2)^2 - norm(r_SCI)^2, 3*r_SCI(2)*r_SCI(3);
    3*r_SCI(1)*r_SCI(3), 3*r_SCI(2)*r_SCI(3), 3*r_SCI(3)^2 - norm(r_SCI)^2];

if flags(2)
    [a_SRP, d_dr_SRP, d_drdot_SRP] = acc_SRP(r_SCI, A, m); % [km/s^2]
else
    a_SRP = zeros(3,1); % [km/s^2]
    d_dr_SRP = zeros(3);
    d_drdot_SRP = zeros(3);
end

if flags(3)
    [a_TB, d_dr_TB, d_drdot_TB] = acc_TB_SCI(r_SCI, et, bodies);
else
    a_TB = zeros(3,1);
    d_dr_TB = zeros(3);
    d_drdot_TB = zeros(3);
end

d_dr = d_dr + d_dr_SRP + d_dr_TB;
d_drdot = d_drdot + d_drdot_SRP + d_drdot_TB;
PHI_dot = [zeros(3),  eye(3);
               d_dr, d_drdot]*PHI;

vec_dot = zeros(6,1);
vec_dot(1:3) = v_SCI;
vec_dot(4:6) = a_SCI + a_SRP + a_TB;
vec_dot(7:42) = reshape(PHI_dot,36,1);
end
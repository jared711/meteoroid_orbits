function [vec_dot] = dynamics_nondim(t, vec, et0, flags, norm_units, A, m, bodies, frame)
%DYNAMICS_NONDIM calculates the time derivative of a non-dimensional
%state vector given in cartesian coordinates in an inertial frame about the
%central body. Can take drag and solar radiation pressure into account.
% 
% [vec_dot] = dynamics_nondim(t, vec, et0, flags, norm_units, A, m, bodies)
% 
% Inputs:   t [NON] (scalar) dimensionless time used by ode78e, ode113 etc...
%           vec [NON;NON;NON] (42x1) dimensionless cartesian state vector
%           et0 [sec] (scalar) ephemeris time (seconds from January 1, 2000)
%           flags [] (4x1) flags = [flag_drag; flag_SRP; flag_TB; flag_J2];
%           norm_units [km, sec] (2x1) normalizing units DU and VU
%           A = [m^2] (scalar) cross-sectional area
%           m = [kg] (scalar) mass of meteoroid
%           bodies [] {1xN} cell array of names of third bodies          
% 
% Outputs:  vec_dot [km/s;km/s^2;NON] (42x1) d/dt of cartesian state vector
% 
% See also: 

% Author: Jared Blanchard 	Date: 2022/01/28 13:58:56 	Revision: 0.1 $

% options (struct)
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

if isrow(vec);    vec = vec';                     end

DU = norm_units(1);
VU = norm_units(2);
TU = DU/VU;

et = et0 + t*TU; % [sec] ephemeris time
% gm = cspice_bodvrd( 'EARTH', 'GM', 1 );
gm = 1; % in dimensionless problems, gm = 1;

x = vec(1:6); % state vector
r = vec(1:3); % position vector
v = vec(4:6); % velocity vector
a = (-gm/(norm(r)^3))*r; % acceleration vector

a_tot = a;
dadr_tot = zeros(3);
dadv_tot = zeros(3);

x_dim = dim(x,DU,VU);
r_dim = x_dim(1:3);

if flags(1) % drag
    if frame == "ECI"
        [a_drag, dadr_drag, dadv_drag] = acc_drag(x_dim, A, m); % [km/s^2]
        a_tot = a_tot + a_drag*DU/VU^2; % *TU^2 is the same as dividing by the acceleration unit (VU^2/DU)
        dadr_tot = dadr_tot + dadr_drag*TU^2;
        dadv_tot = dadv_tot + dadv_drag*TU;
    end
end

if flags(2) % SRP
    if frame == "ECI"
        r_earth_ECI = cspice_spkpos('EARTH', et, 'J2000', 'NONE', 'SUN');
        [a_SRP, dadr_SRP, dadv_SRP] = acc_SRP(r_dim + r_earth_ECI, A, m); % [km/s^2]
    else
        [a_SRP, dadr_SRP, dadv_SRP] = acc_SRP(r_dim, A, m); % [km/s^2]
    end
    a_tot = a_tot + a_SRP*DU/VU^2;
    dadr_tot = dadr_tot + dadr_SRP*TU^2;
    dadv_tot = dadv_tot + dadv_SRP*TU;
end

if flags(3) % Third Body
    if frame == "ECI"
        [a_TB, dadr_TB, dadv_TB] = acc_TB_ECI(r_dim, et, bodies);
    elseif frame == "SCI"
        [a_TB, dadr_TB, dadv_TB] = acc_TB_SCI(r_dim, et, bodies);
    end
    a_tot = a_tot + a_TB*DU/VU^2;
    dadr_tot = dadr_tot + dadr_TB*TU^2;
    dadv_tot = dadv_tot + dadv_TB*TU;
end

try 
    if flags(4) % J2
        if frame == "ECI"
            [a_J2, dadr_J2, dadv_J2] = acc_J2(x_dim);
            a_tot = a_tot + a_J2*DU/VU^2;
            dadr_tot = dadr_tot + dadr_J2*TU^2;
            dadv_tot = dadv_tot + dadv_J2*TU;
        end
    end
catch
end
        
vec_dot = zeros(6,1);
vec_dot(1:3) = v;
vec_dot(4:6) = a_tot;

if length(vec) == 42
    PHI = reshape(vec(7:42),6,6);

    % partials from central gravity field (so I don't have to compute this
    % when there is no STM
    dadr_tot = dadr_tot + gm/(norm(r)^5)*...
       [3*r(1)^2 - norm(r)^2,          3*r(1)*r(2),          3*r(1)*r(3);
                 3*r(1)*r(2), 3*r(2)^2 - norm(r)^2,          3*r(2)*r(3);
                 3*r(1)*r(3),          3*r(2)*r(3), 3*r(3)^2 - norm(r)^2];
    
    % variational equations
    PHI_dot = [zeros(3),   eye(3);
               dadr_tot, dadv_tot]*PHI;
    vec_dot(7:42) = reshape(PHI_dot,36,1);
elseif length(vec) == 6
else
    warning("input vector should be 6x1 or 42x1")
end

end
function [a_drag, d_dr, d_drdot] = acc_drag(x_ECI, A, m)
%ACC_DRAG computes acceleration due to drag in the ECI frame
% 
% [a_drag, drddot] = ACC_DRAG(x_ECI, A, m)
% 
% Inputs:   x_ECI [km;km/s] (6x1) state in ECI frame
%           A = [m^2] (scalar) cross-sectional area
%           m = [kg] (scalar) mass of meteoroid
% 
% Outputs:  a_drag [km/s^2] (3x1) acceleration due to drag
%           d_dr [] (3x3) partial derivatives w.r.t r
%           d_drdot [] (3x3) partial derivatives w.r.t rdot
% 
% See also: 

% Author: Jared Blanchard 	Date: 2022/02/03 09:05:21 	Revision: 0.1 $

R_E = cspice_bodvrd( 'EARTH', 'RADII', 3); %[km] radii of Earth Ellipsoid
r_ECI = x_ECI(1:3);
v_ECI = x_ECI(4:6);
h = abs(norm(r_ECI) - R_E(1)); %[km]
if h <= 180 % [km] vida, D., Gural, P. S., Brown, P. G., Campbell-Brown, M., & Wiegert, P. (2020). Estimating trajectories of meteors: An observational Monte Carlo approach - I. Theory. Monthly Notices of the Royal Astronomical Society, 491(2), 2688–2705. https://doi.org/10.1093/mnras/stz3160
    % calculate atmospheric drag assuming the atmosphere is moving along
    % with the surface of the central body
    omega = 7.292115e-5; %[rad/s]
    H = 8.5e3; %[m]
    c_d = 0.5; % coefficient of drag (2/3/22 may need to vary with altitude)
%     c_d = 1;
    ang_vel = [0;0;omega]; % [rad/s] angular velocity vector of central body  x
    v_atm = cross(ang_vel,r_ECI); % [km/s] velocity vector of atmosphere  x
    v_rel = (v_ECI - v_atm)*1000; % [m/s] velocity of meteor relative to atmosphere x 
%     m = options.m/1000; % [kg] mass of meteoroid x
%     c_d = options.c_d; % 2.3; % [] coefficient of drag of meteoroid x
%     A = options.A; % [m^2] x
%     H = options.H; % [m] scale height x
    rho_0 = 1.255; % [kg/m^3] density of atmosphere at sea level x 

    rho = rho_0*exp(-((h*1000)/H)); %[kg/m^3] % density of atmosphere at altitude %% ToDo change this to GD altitude
%         set_rho_atm
%         rho = interp1(rho_x,rho_1.rho,(r_mag-R_E(1)*1000)/1000); %[kg/m^3]
% 2/3/22 I may need to change the drag coefficient from 0.5 to 1 as it goes up in the atmosphere 


    a_drag = -1/2*c_d*(A/m)*rho*norm(v_rel)*v_rel; % [m/s^2]
    a_drag = a_drag/1000; % [km/s^2]
%     a_drag = R_SCI2ECI(1:3,1:3)\a_drag; % convert back to SCI frame (inverse of SCI2ECI)

%         a_drag = -2*c_d*rho*norm(v_rel)^2*(1.21/m)*(m/meteor.pd3d)^(2/3)*unit(v_rel)/1000; % Vida, D., Brown, P. G., & Campbell-Brown, M. (2018). Modelling the measurement accuracy of pre-atmosphere velocities of meteoroids. Monthly Notices of the Royal Astronomical Society, 479(4), 4307–4319. https://doi.org/10.1093/mnras/sty1841
    d_drdot = -1/2*c_d*(A/m)*rho*((v_rel*v_rel')/norm(v_rel) + norm(v_rel)*eye(3));
    d_dr = -1/2*c_d*(A/m)*norm(v_rel)*v_rel * (-1000*rho_0)/H*exp(-h*1000/H) - d_drdot * [0, -omega, 0; omega, 0, 0; 0, 0, 0];
else
    a_drag = zeros(3,1); % [km/s^2]
    d_dr = zeros(3);
    d_drdot = zeros(3);
end



end

function [a_TB, d_dr, d_drdot] = acc_TB_ECI(r_ECI, et, bodies)
%ACC_TB_ECI computes acceleration due to third body effects in the ECI frame
% 
% [a_TB, d_dr, d_drdot] = ACC_TB_ECI(r_ECI, et, bodies)
% 
% Inputs:   r_ECI [km] (3x1) position in ECI frame
%           et = [sec] (scalar) ephemeris time
%           bodies = [] {Nx1} cell array containing names of bodies to
%           include
% 
% Outputs:  a_TB [km/s^2] (3x1) acceleration due to third-body effects
%           d_dr [] (3x3) partial derivatives w.r.t r
%           d_drdot [] (3x3) partial derivatives w.r.t rdot
% 
% See also: acc_J2, acc_SRP, acc_drag, acc_TB_SCI

% Author: Jared Blanchard 	Date: 2022/02/03 09:28:32 	Revision: 0.1 $

a_TB = zeros(3,1);
d_dr = zeros(3);
d_drdot = zeros(3);
for body = bodies
    r_body = cspice_spkpos(body{1}, et, 'J2000', 'NONE', 'EARTH');
    r_i = r_ECI - r_body;
    gm = cspice_bodvrd(body{1}, 'GM', 1 );
    a_TB = a_TB + (-gm/(norm(r_i)^3))*r_i;
    d_dr = d_dr - gm*(eye(3)/norm(r_i)^3 - 3*(r_i*r_i')/norm(r_i)^5);
end

end

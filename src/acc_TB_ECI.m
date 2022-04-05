function [a_TB, d_dr, d_drdot] = acc_TB_ECI(r_ECI, et, bodies)
%ACC_TB Summary of this function goes here
% 
% [a_TB, d_dr, d_drdot] = ACC_TB_ECI(r_ECI, et, options)
% 
% Inputs: 
% 
% Outputs: 
% 
% See also: 

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

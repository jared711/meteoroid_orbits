function [x_ECEF] = latlonazel2ECEF(meteor)
%LATLONAZEL2ECEF Summary of this function goes here
% 
% [OUTPUTARGS] = LATLONAZEL2ECEF(INPUTARGS)
% 
% Inputs: 
% 
% Outputs: 
% 
% See also: 

% Author: Jared Blanchard 	Date: 2022/01/28 17:28:35 	Revision: 0.1 $

radii = cspice_bodvrd( 'EARTH', 'RADII', 3 );
flat = (radii(1) - radii(3))/radii(1);
r_ECEF = cspice_georec(deg2rad(meteor.lon), deg2rad(meteor.lat), meteor.h, radii(1), flat);


v_ENU = azel2cart(meteor.az, meteor.el, -meteor.v);
R = rotzd(90 + meteor.lon)*rotxd(90 - meteor.lat);
v_ECEF = R*v_ENU;


x_ECEF = [r_ECEF; v_ECEF];

end

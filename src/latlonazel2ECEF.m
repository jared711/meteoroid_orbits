function [x_ECEF, J2] = latlonazel2ECEF(meteor, rand)
%LATLONAZEL2ECEF Summary of this function goes here
% 
% [x_ECEF, J2] = LATLONAZEL2ECEF(meteor, rand)
% 
% Inputs: 
% 
% Outputs: 
% 
% See also: 

% Author: Jared Blanchard 	Date: 2022/01/28 17:28:35 	Revision: 0.1 $
if nargin < 2;  rand = 0;   end

if isstruct(meteor)
    if rand
        lat = meteor.lat + meteor.lat_err*randn;
        lon = meteor.lon + meteor.lon_err*randn;
        h = meteor.h + meteor.h_err*randn;
        az = meteor.az + meteor.az_err*randn;
        el = meteor.el + meteor.el_err*randn;
        v = meteor.v + meteor.v_err*randn;
    else
        lat = meteor.lat;
        lon = meteor.lon;
        h = meteor.h;
        az = meteor.az;
        el = meteor.el;
        v = meteor.v;
    end
else
    lat = meteor(1);
    lon = meteor(2);
    h = meteor(3);
    az = meteor(4);
    el = meteor(5);
    v = meteor(6);
end
    
radii = cspice_bodvrd( 'EARTH', 'RADII', 3 );
flat = (radii(1) - radii(3))/radii(1);
r_ECEF = cspice_georec(deg2rad(lon), deg2rad(lat), h, radii(1), flat);

v_ENU = azel2cart(az, el, -v);
R = rotzd(90 + lon)*rotxd(90 - lat);
v_ECEF = R*v_ENU;

x_ECEF = [r_ECEF; v_ECEF];

J2 = eye(6); % 2/7/22 To-Do implement this Jacobian

end

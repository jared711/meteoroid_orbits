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
% r_ECEF = GD2ECEF(lat,lon,h); % Replace custom function with Spice equivalent 

v_ENU = azel2cart(az, el, -v);
R = rotzd(90 + lon)*rotxd(90 - lat);
v_ECEF = R*v_ENU;

x_ECEF = [r_ECEF; v_ECEF];

r_E = 6378.137; % [km]
e_E = 0.081819190842621; % []

N = r_E/sqrt(1-e_E^2*sind(lat)^2);
dNdlat = r_E*e_E^2*sind(lat)*cosd(lat)/sqrt(1-e_E^2*sind(lat)^2)^3;

S1 = [-(N+h)*sind(lat)*cosd(lon) + dNdlat*cosd(lat)*cosd(lon),  -(N+h)*cosd(lat)*sind(lon), cosd(lat)*cosd(lon);
      -(N+h)*sind(lat)*sind(lon) + dNdlat*cosd(lat)*sind(lon),   (N+h)*cosd(lat)*cosd(lon), cosd(lat)*sind(lon);
     (N*(1-e_E^2) + h)*cosd(lat) + dNdlat*(1-e_E^2)*sind(lat),                           0,           sind(lat)];
S2 = zeros(3);
dRdlat = rotzd(lon+90)*[0,             0,            0;
                        0,  sind(90-lat), cosd(90-lat);
                        0, -cosd(90-lat), sind(90-lat)];
dRdlon = [-sind(90+lon), -cosd(90+lon), 0;
           cosd(90+lon), -sind(90+lon), 0;
                      0,             0, 0]*rotxd(90-lat);
S31 = dRdlat*[v*cosd(el)*sind(az);
              v*cosd(el)*cosd(az);
                       v*sind(el)];
S32 = dRdlon*[v*cosd(el)*sind(az);
              v*cosd(el)*cosd(az);
                       v*sind(el)];
S33 = zeros(3,1);
S3 = [S31, S32, S33];
S4 = R*[v*cosd(el)*cosd(az), -v*sind(el)*sind(az), cosd(el)*sind(az);
       -v*cosd(el)*sind(az), -v*sind(el)*cosd(az), cosd(el)*cosd(az);
                          0,           v*cosd(el),          sind(el)];
J2 = [S1, S2;
      S3  S4];
                           
end

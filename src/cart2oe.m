function [a,e,i,OMEGA,omega,nu,ang,ME] = cart2oe(x, mu)
%CART2OE converts cartesian coordinates in inertial frame to orbital
%elements. This function is able to handle any orbit type , including 
% equatorial and circular orbits. The output extra_angle will output the
% following values for each orbit type
% 1) Elliptical equatorial
% Longitude of periapsis , ang = Pi = OMEGA + omega
% 2) Circular inclined
% Argument of latitude , ang = u = omega + nu
% 3) Circular equatorial
% True latitude , ang = l = OMEGA + omega + nu
% Otherwise , ang will be undefined ( NaN )
% 
% [a, e, i, OMEGA, omega, nu, ang, ME] = CART2OE(rv_ECI, mu)
% 
% Inputs:   x = [km; km/s] (6x1) state in body-centered inertial frame
%           mu = [km^3/s^2] (scalar) gravitational parameter of body
% 
% Outputs:  a = [km] (scalar) semimajor axis
%           e = [] (scalar) eccentricity
%           i = [deg/rad] (scalar) inclination
%           OMEGA = [deg/rad] (scalar) right ascension of the ascending node
%           omega = [deg/rad] (scalar) argument of periapsis
%           nu = [deg/rad] (scalar) true anomaly
%           ang [deg/rad] (scalar) extra angle for special cases
%           ME = [km^2/s^2] (scalar) mechanical energy
% 
% See also: 

% Author: Jared Blanchard 	Date: 2021/02/15 09:31:37 	Revision: 0.1 $
% Adapted from ECI2OE.m used in Stanford's AA279A course.
% Author: Jared Blanchard 	Date: 2022/01/28 012:25:10 	Revision: 0.2 $

if isrow(x);        x = x';                         end
if ~iscolumn(x);    error('x should be a vector');  end

rVec = x(1:3);
vVec = x(4:6);
r = norm(rVec);
v = norm(vVec);
% Create all necessary vectors
hVec = cross(rVec, vVec);
h = norm(hVec);
nVec = cross([0, 0, 1], hVec);
n = norm(nVec);
eVec = (1/mu)*((v^2 - mu/r)*rVec - dot(rVec, vVec)*vVec);
e = norm(eVec);
% Compute the size of the orbit
ME = 0.5* v^2 - mu/r;
if e ~= 1
    a = -mu/(2*ME);
    p = a*(1 - e ^2);
else
    a = Inf ;
    p = h^2/mu;
end
% Compute the orientation of the orbit
i = acosd (hVec(3)/h);
OMEGA = acosd (nVec(1)/n);
omega = acosd(dot(nVec, eVec)/(n*e));
nu = acosd(dot(eVec, rVec)/(e*r));
% Place angles in the correct domains
if nVec (2) < 0
    OMEGA = 360 - OMEGA;
end
if eVec (3) < 0
    omega = 360 - omega;
end
if dot(rVec, vVec) < 0
    nu = 360 - nu ;
end
% Account for any special cases
if i == 0 && e ~= 0 % Elliptical equatorial
    % Provide the longitude of periapsis (PI = Om + w)
    ang = acosd ( eVec (1) /e);
if eVec (2) < 0
    ang = 360 - ang;
end
elseif i ~= 0 && e == 0 % Circular inclined
    % Provide the argument of latitude (u = w + anom )
    ang = acosd ( dot (nVec , rVec )/(n*r));
if rVec (3) < 0
    ang = 360 - ang;
end
elseif i == 0 && e == 0 % Circular equatorial
    % Provide the true latitude ( l = Om + w + anom )
    ang = acosd ( rVec (1) /r);
if rVec (2) < 0
    ang = 360 - ang;
end
else
    % Default output for ang
    ang = NaN ;
end

end
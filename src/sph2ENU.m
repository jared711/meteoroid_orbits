function [x_ENU, J1] = sph2ENU(x_sph)
%SPH2ENU Converts a 6D state vector from spherical coordinates
%[r;az;el;rdot,azdot,eldot] to cartesian coordinates and outputs the
%jacobian matrix
% 
% [x_ENU, J] = SPH2ENU(x_sph)
% 
% Inputs:   x_sph {6x1} [km;deg;deg;km/s;deg/s;deg/s] spherical state 
%               vector of form [r, az, el, rdot, azdot, eldot]'
% 
% Outputs:  x_ENU {6x1} [km;km/s] cartesian state vector
%           J {6x6} [] Jacobian computed at x_sph
% 
% Author: Jared Blanchard 	Date: 2021/06/03 19:02:34 	Revision: 0.1 $

r     = x_sph(1);   az     = x_sph(2);  el     = x_sph(3);
r_dot = x_sph(4);   az_dot = x_sph(5);  el_dot = x_sph(6);

r_ENU = r*[cosd(el)*sind(az);
           cosd(el)*cosd(az);
           sind(el)];
     
v_ENU = [r_dot*cosd(el)*sind(az) + r*deg2rad(az_dot)*cosd(el)*cosd(az) - r*deg2rad(el_dot)*sind(el)*sind(az);
         r_dot*cosd(el)*cosd(az) - r*deg2rad(az_dot)*cosd(el)*sind(az) - r*deg2rad(el_dot)*sind(el)*cosd(az);
         r_dot*sind(el) + r*deg2rad(el_dot)*cosd(el)];
     
x_ENU = [r_ENU; v_ENU];

J11 = [cosd(el)*sind(az),  r*cosd(el)*cosd(az), -r*sind(el)*sind(az);
       cosd(el)*cosd(az), -r*cosd(el)*sind(az), -r*sind(el)*cosd(az);
                sind(el),                    0,           r*cosd(el)];
            

J21 = [deg2rad(az_dot)*cosd(el)*cosd(az)  - deg2rad(el_dot)*sind(el)*sind(az),...
       r_dot*cosd(el)*cosd(az)   - r*deg2rad(az_dot)*cosd(el)*sind(az) - r*deg2rad(el_dot)*sind(el)*cosd(az),...
       -r_dot*sind(el)*sind(az)  - r*deg2rad(az_dot)*sind(el)*cosd(az) - r*deg2rad(el_dot)*cosd(el)*sind(az);
       -deg2rad(az_dot)*cosd(el)*sind(az) - deg2rad(el_dot)*sind(el)*cosd(az),...
       r_dot*cosd(el)*sind(az)   - r*deg2rad(az_dot)*cosd(el)*cosd(az) + r*deg2rad(el_dot)*sind(el)*sind(az),...
       -r_dot*sind(el)*cosd(az)  + r*deg2rad(az_dot)*sind(el)*sind(az) - r*deg2rad(el_dot)*cosd(el)*cosd(az);
       deg2rad(el_dot)*cosd(el),...
       0,...
       r_dot*cosd(el) - r*deg2rad(el_dot)*sind(el)];
           
J1 = [J11, zeros(3);
      J21,     J11];
              
end

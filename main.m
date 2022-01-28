% This should be the only script/function that people need to deal with

addpath(genpath('src'))
cspice_furnsh('src/kernels/standard.tm') % load the SPICE kernels

% x_sph = ones(6,1);
% Sigma0 = eye(6);
% 
% [x_ENU, J1] = sph2ENU(x_sph);
% Sigma1 = J1*Sigma0*J1';

x_ENU = [28.3109;8.3590;109.7792;10.1808;24.5716;64.8084];
time = [2008,1,2,18,10,20.5249];
Sigma1 = eye(6);

% geodetic location of ALTAIR
phi = 0.1640; %[rad] latitude
lambda = 2.9231; %[rad] longitude
h = 0.0835; % [km] altitude

[x_SCI, S] = ENU2SCI(x_ENU, phi, lambda, h, time);
Sigma2 = S*Sigma1*S';



GM = 1.327124400179870e11; %[km^3/s^2]
AU = 149597927; %[km]
[a,e,i,Omega,omega] = cart2oe(x_SCI, GM);
a = a/AU

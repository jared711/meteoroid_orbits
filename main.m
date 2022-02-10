% This should be the only script/function that people need to deal with
clc
clear all
% close all

addpath(genpath('src'))
cspice_furnsh('src/kernels/loadkernels.tm') % load the SPICE kernels

flag_drag = true;
flag_SRP = true;
flag_TB = true;
flags = [flag_drag; flag_SRP; flag_TB];
bodies = {'EARTH','MOON','JUPITER_BARYCENTER'};
tf = -86400*10*365; % [sec] 86400 seconds in a day
evfcn_e = @(t,x) ef_rval3d(t,x,1e6); % past Earth's SOI
evfcn_j = @(t,x) ef_rval3d(t,x,1e9); % past Jupiter's Orbit
tol = 1e-12;
trace = 0;

% % options.TB = {'EARTH','MOON','JUPITER_BARYCENTER'};
% % options.time0 = time;
% % % options.SRP = 'true';
% % options.c_srp = 1;
% % % options.m = 1e-6; %[g]
% % options.m = kosice.m; %[g]
% % % options.A = 4.1347e-09; % [m^2] % effective area for SRP and drag
% % options.A = kosice.A; % [m^2]
% % options.drag = 'true';

% x_sph = ones(6,1);
% Sigma0 = eye(6);
% 
% [x_ENU, J1] = sph2ENU(x_sph);
% Sigma1 = J1*Sigma0*J1';

%%%%%%%%%%%%%% these are the same as the following meteors except I changed
%%%%%%%%%%%%%% a line of code in analyzeheads07.m line 144 to not have a
%%%%%%%%%%%%%% negative sign
% x_ENU = [28.3109;8.3590;109.7792;-23.1549;14.7290;-64.4550];
% time = [2008,1,2,18,10,20.5249];
% x_ENU = [30.1501;7.4197;109.0860;-12.6312;-5.4119;-57.3434];
% time = [2008,1,2,18,10,22.5422];
%%%%%%%%%%%%%%%%

load altairMeteors.mat
rand = 0;
num = 8;
N = 1;
oes = zeros(N,5);
for j = 1:N
    if rand
        x_sph = [meteors{num}.range + meteors{num}.err_range*randn;
                meteors{num}.fit_az + meteors{num}.bore_az + meteors{num}.err_az*randn;
                meteors{num}.fit_el + meteors{num}.bore_el + meteors{num}.err_el*randn;
                meteors{num}.range_rate + meteors{num}.err_range_rate*randn;
                meteors{num}.fit_azdot + meteors{num}.err_azdot*randn;
                meteors{num}.fit_eldot + meteors{num}.err_eldot*randn];
    else
        x_sph = [meteors{num}.range;
                 meteors{num}.fit_az + meteors{num}.bore_az;
                 meteors{num}.fit_el + meteors{num}.bore_el;
                 meteors{num}.range_rate;
                 meteors{num}.fit_azdot;
                 meteors{num}.fit_eldot];
    end
    x_ENU = sph2ENU(x_sph);
    time = meteors{num}.time;


    % % x_ENU = [28.3109;8.3590;109.7792;10.1808;24.5716;64.8084];
    % time = [2008,1,2,18,10,20.5249];
    % x_ENU = [30.1501;7.4197;109.860;18.6647;2.2897;55.8882];
    % time = [2008,1,2,18,10,22.5422];
    et0 = cspice_str2et(datestr(time));
    m = 1e-9; % [kg]
    A = 4.1347e-9; % [m^2]

    % Sigma0 = eye(6);

    % geodetic location of ALTAIR
    phi = 0.1640; %[rad] latitude
    lambda = 2.9231; %[rad] longitude
    h = 0.0835; % [km] altitude
    [x_ECEF, R1] = ENU2ECEF(x_ENU, phi, lambda, h);
    % Sigma1 = R1*Sigma0*R1';

    % [x_SCI, S] = ENU2SCI(x_ENU, phi, lambda, h, time);
    % Sigma2 = S*Sigma1*S';

    % verification_data
    % x_ECEF = latlonazel2ECEF(kosice);
    % A = kosice.A;
    % m = kosice.m;
    % time = kosice.time;
    % et0 = cspice_str2et(datestr(time));

    [x_ECI, R2] = ECEF2ECI(x_ECEF, et0);
    % Sigma2 = R2*Sigma1*R2';

    vec0_ECI = [x_ECI; reshape(eye(6),36,1)];
    [tt_ECI, xx_ECI] = ode78e(@(t,vec) dynamics_ECI(t, vec, et0, flags, A, m, {'MOON','SUN'}), 0, tf, vec0_ECI, tol, trace, evfcn_e);
    x_ECI_f = xx_ECI(end,1:6)';
    % R3 = reshape(xx_ECI(end,7:42),6,6);
    t_int_ECI = tt_ECI(end);
    % Sigma3 = R3*Sigma2*R3';

    [x_SCI, R4] = ECI2SCI(x_ECI_f, et0 + t_int_ECI);
    % Sigma4 = R4*Sigma3*R4';

    vec0_SCI = [x_SCI; reshape(eye(6),36,1)];
    [tt_SCI, xx_SCI] = ode78e(@(t,x) dynamics_SCI(t, x, et0 + t_int_ECI, flags, A, m, {'MOON','EARTH','JUPITER_BARYCENTER'}), 0, tf, vec0_SCI, tol, trace, evfcn_j);
    x_SCI_f = xx_SCI(end,1:6)';
    % R5 = reshape(xx_SCI(end,7:42),6,6);
    t_int_SCI = tt_SCI(end);
    % Sigma5 = R5*Sigma4*R5';

    GM = cspice_bodvrd( 'SUN', 'GM', 1 );
    AU = cspice_convrt(1, 'AU', 'km'); %[km]

    [a,e,i,Omega,omega] = cart2oe(x_SCI_f, GM);
    a = a/AU;
%     sprintf('%1.4f & %1.4f & %1.4f & %1.4f & %1.4f',a,e,i,Omega,omega)
    oes(j,:) = [a,e,i,Omega,omega];
    
end

% std(oes)

% kosice.oe
% let's plot how the variables change over time

% figure('WindowState','maximized')
figure('Position',[1228         701         829         695])
plot(xx_SCI(:,1)/AU,xx_SCI(:,2)/AU,'k','LineWidth',2)
hold on
plot_circle(1,'b',':',4,'on')
plot_circle(5.2,'r',':',4,'on')
xlabel('X [AU]')
ylabel('Y [AU]')
legend('Meteor Orbit','Earth Orbit','Jupiter Orbit')
% savefigs('altair2all')
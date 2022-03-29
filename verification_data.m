% Dmitriev, V., Lupovka, V., & Gritsevich, M. (2015). Orbit determination based on meteor observations using numerical integration of equations of motion. Planetary and Space Science, 223–235. https://doi.org/10.1016/j.pss.2015.06.015

% Oijarvi
oijarvi.lat = 64.78;    oijarvi.lat_err = 0.02; %[deg]
oijarvi.lon = 26.91;    oijarvi.lon_err = 0.01; %[deg]
oijarvi.h = 77.00;      oijarvi.h_err = 1; %[km]
oijarvi.az = 156.2;     oijarvi.az_err = 0.3; %[deg]
oijarvi.el = 25.80;     oijarvi.el_err = 0.3; %[deg]
oijarvi.v = 13.8;       oijarvi.v_err = 1; %[km/s]
oijarvi.A = [];
oijarvi.m = [];
oijarvi.time = [2010,12,26,14,06,09];   oijarvi.time_err = 5; %[sec]
oijarvi.name = 'Oijarvi';
oijarvi.oe = [2.4582, 0.6011, 2.8018, 94.5264, 352.7669];
oijarvi.Sigma = diag([oijarvi.lat_err^2, oijarvi.lon_err^2, oijarvi.h_err^2, oijarvi.az_err^2, oijarvi.el_err^2, oijarvi.v_err^2]);

% Mikkeli
mikkeli.lat = 61.46;    mikkeli.lat_err = 0.01; %[deg]
mikkeli.lon = 26.90;    mikkeli.lon_err = 0.01; %[deg]
mikkeli.h = 82.10;      mikkeli.h_err = 0.5; %[km]
mikkeli.az = 238.94;    mikkeli.az_err = 0.4; %[deg]
mikkeli.el = 55.06;     mikkeli.el_err = 0.2; %[deg]
mikkeli.v = 14.98;      mikkeli.v_err = 0.1; %[km/s]
mikkeli.A = [];
mikkeli.m = [];
mikkeli.time = [2013,09,13,22,33,37];   mikkeli.time_err = 5; %[sec]
mikkeli.name = 'Mikkeli';
mikkeli.oe = [1.4354, 0.3652, 12.2108, 171.1334, 229.6174];
mikkeli.Sigma = diag([mikkeli.lat_err^2, mikkeli.lon_err^2, mikkeli.h_err^2, mikkeli.az_err^2, mikkeli.el_err^2, mikkeli.v_err^2]);

% Annama
annama.lat = 67.93;    annama.lat_err = 0.01; %[deg]
annama.lon = 30.76;    annama.lon_err = 0.01; %[deg]
annama.h = 83.90;      annama.h_err = 0.5; %[km]
annama.az = 176.1;     annama.az_err = 0.5; %[deg]
annama.el = 34.32;     annama.el_err = 0.5; %[deg]
annama.v = 24.21;      annama.v_err = 0.5; %[km/s]
annama.A = [];
annama.m = [];
annama.time = [2014,04,18,22,14,09.3];   annama.time_err = 0.5; %[sec]
annama.name = "Annama";
annama.oe = [2.0007, 0.6825, 14.6109, 28.6098, 264.5874];
annama.Sigma = diag([annama.lat_err^2, annama.lon_err^2, annama.h_err^2, annama.az_err^2, annama.el_err^2, annama.v_err^2]);

% Annama
haapavesi.lat = 66.52;   haapavesi.lat_err = 0.03; %[deg]
haapavesi.lon = 25.16;    haapavesi.lon_err = 0.02; %[deg]
haapavesi.h = 70.95;      haapavesi.h_err = 1; %[km]
haapavesi.az = 357.25;     haapavesi.az_err = 0.3; %[deg]
haapavesi.el = 11.05;     haapavesi.el_err = 0.3; %[deg]
haapavesi.v = 14.78;      haapavesi.v_err = 0.3; %[km/s]
haapavesi.A = [];
haapavesi.m = [];
haapavesi.time = [2014,09,25,03,12,15.0];   haapavesi.time_err = 1; %[sec]
haapavesi.name = "Haapavesi";
haapavesi.oe = [2.5301, 0.6042, 9.2407, 181.8271, 174.6757];
haapavesi.Sigma = diag([haapavesi.lat_err^2, haapavesi.lon_err^2, haapavesi.h_err^2, haapavesi.az_err^2, haapavesi.el_err^2, haapavesi.v_err^2]);

% Buzzard Coulee
buzzard.lat = 53.169;   buzzard.lat_err = 0.001; %[deg]
buzzard.lon = -10.059;    buzzard.lon_err = 0.001; %[deg]
buzzard.h = 63.8;      buzzard.h_err = 0.7; %[km]
buzzard.az = 347.5;     buzzard.az_err = 0.4; %[deg]
buzzard.el = 66.7;     buzzard.el_err = 0.4; %[deg]
buzzard.v = 18.0;      buzzard.v_err = 0.4; %[km/s]
buzzard.A = [];
buzzard.m = [];
buzzard.time = [2008,11,21,00,26,43.0];   buzzard.time_err = 1; %[sec]
buzzard.name = "Buzzard Coulee";
buzzard.oe = [1.25, 0.23, 25.0, 238.93739, 211.3];
buzzard.Sigma = diag([buzzard.lat_err^2, buzzard.lon_err^2, buzzard.h_err^2, buzzard.az_err^2, buzzard.el_err^2, buzzard.v_err^2]);

% Jansen-Sturgeon, T., Sansom, E. K., & Bland, P. A. (2019). Comparing analytical and numerical approaches to meteoroid orbit determination using Hayabusa telemetry. Meteoritics and Planetary Science, 54(9), 2149–2162. https://doi.org/10.1111/maps.13376
hayabusa1.lat = -29.0243;   hayabusa1.lat_err = 0;% [deg]   
hayabusa1.lon = 131.1056;   hayabusa1.lon_err = 0; % [deg] 
hayabusa1.h = 99.880;       hayabusa1.h_err = 0; % [km]
hayabusa1.az = 290.5220;    hayabusa1.az_err = 0; % [deg]
hayabusa1.el = 10.0173;     hayabusa1.el_err = 0; % [deg]     
hayabusa1.v = 11.7251;      hayabusa1.v_err = 0.01; % [km/s]       
hayabusa1.time = [2010,06,13,13,51,56.6];
hayabusa1.name = 'Hayabusa Spacecraft';
hayabusa1.m = 415; % [kg]
hayabusa1.A = 2.15; % [m^2]
hayabusa1.oe = [1.32381, 0.25732, 1.68383, 82.46569, 147.47773];


hayabusa2.lat = -29.6545; % [deg]
hayabusa2.lon = 133.0768; % [deg] 
hayabusa2.h = 64.710; % [km]
hayabusa2.az = 289.2733; % [deg]
hayabusa2.el = 8.7955; % [deg]     
hayabusa2.v = 11.3305; % [km/s]       
hayabusa2.time = [2010,6,13,13,52,16.0];
hayabusa2.name = 'hayabusa capsule';
hayabusa2.m = 20; % [kg]
hayabusa2.A = 0.126; % [m^2]
hayabusa2.oe = [1.32381, 0.25732, 1.68383, 82.46569, 147.47773];

kosice.m = 3500; %[kg]
kosice.rho = 3.4; %[g/cm^3]
kosice.r = (3*kosice.m*1000/(4*pi*kosice.rho))^(1/3); % [cm]
kosice.A = pi*kosice.r^2; %[cm^2]
kosice.A = kosice.A/1e4; %[m^2]
kosice.lat = 48.467;    kosice.lat_err = 0.021; %[deg]
kosice.lon = 20.705;    kosice.lon_err = 0.011;%[deg]
kosice.h = 68.3;        kosice.h_err = 1.4;%[km]
kosice.az = 252.6;      kosice.az_err = 4; %[deg]
kosice.el = 58.8;       kosice.el_err = 2; %[deg]
kosice.v = 15.0;        kosice.v_err = 0.3; %[km/s]
kosice.time = [2010,02,28,22,24,47.0];
kosice.name = 'Kosice';
kosice.oe = [2.725, 0.649, 2.015, 340.14579, 204.108];
kosice.Sigma = diag([kosice.lat_err^2, kosice.lon_err^2, kosice.h_err^2, kosice.az_err^2, kosice.el_err^2, kosice.v_err^2]);

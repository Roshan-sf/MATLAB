%% Problem 8.16 (Non-Hohmann Cruise)
clc
clear

muSun = 1.32712*10^11; %km3/s2
%muEarth = 398600; %km3/s2
%muMars = 42828; %km3/s2
rS2E = 149.6*10^6; %km
rS2M = 227.9*10^6; %km

r1 = 190 + 6378; % km
r2p = 300 + 3390; % km
JD2000 = 2451545.0; 
JDstart = juliandate(2005, 8, 15, 12, 0, 0);
JDend = juliandate(2006, 3, 15, 12, 0, 0);

% Where are the planets?
T0e = (JDstart - JD2000)/36525;
T0m = (JDend - JD2000)/36525;
disp("T0 is positive as the start date is after the year 2000")

[Earth_coes] = AERO351planetary_elements2(3, T0e);
[Mars_coes] = AERO351planetary_elements2(4, T0m);

% Earth
omegaE = Earth_coes(5) - Earth_coes(4); % deg
MeE = Earth_coes(6) - Earth_coes(5); % deg
thetaE = Me2TA(deg2rad(MeE), Earth_coes(2)); % rad
h1 = sqrt(rS2E*muSun*(1+Earth_coes(2)*cos(thetaE)));
% h shouldn't be that high for LEO, maybe in GEO??? cuz parking orbit

[rECIe, vECIe] = coes2rv(h1, Earth_coes(2), deg2rad(Earth_coes(3)), ...
    deg2rad(Earth_coes(4)), deg2rad(omegaE), ...
    thetaE, muSun);

disp("Heart: norm of r and v vectors (EARTH) " + ...
    "in PERI matches norm in ECI frame")

% Mars
omegaM = Mars_coes(5) - Mars_coes(4); % deg
MeM = Mars_coes(6) - Mars_coes(5); % deg
thetaM = Me2TA(deg2rad(MeM), Mars_coes(2)); % rad
h2 = sqrt(rS2M*muSun*(1+Mars_coes(2)*cos(thetaM)));

[rECIm, vECIm] = coes2rv(h2, Mars_coes(2), deg2rad(Mars_coes(3)), ...
    deg2rad(Mars_coes(4)), deg2rad(omegaM), ...
    thetaM, muSun);

disp("Heart: norm of r and v vectors (MARS) " + ...
    "in PERI matches norm in ECI frame")

% Lamberts
dt = (JDend - JDstart) * 24 * 60 * 60; %sec
[v1, v2] = lamberts(rECIe, rECIm, dt, muSun);
dv1 = norm(v1 - vECIe);
dv2 = norm(v2 - vECIm);
dvtot = dv1+dv2;

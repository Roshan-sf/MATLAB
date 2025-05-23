%% Roshan Jaiswal-Ferri
%Section - 01 
%Aero 351 Final Exam Question 2: 12/07/24

%% Workspace Prep

format long     %Allows for more accurate decimals
close all;      %Clears all
clear all;      %Clears Workspace
clc;            %Clears Command Window

%% Constant/Global Vars

options = odeset('RelTol',1e-8,'AbsTol',1e-8);
muSun = 1.327e11; %mu values from curtis
muMars = 42828; 
mu = 398600; %earth
muJup = 126686534;
Rmars = 3396; %km
Rearth = 6378; %km
Rjup = 71490;
tol = 1e-8;
ap = 149598000; %1 AU in km
P = (0.75*365.25)*86400; %period in seconds
Re = 149598000;
flyr = Rearth + 10000; %km

%% Finding COEs of Elliptical Orbit

a = ((muSun*P^2)/(4*pi^2))^(1/3);
ecc = (ap-a)/a;
theta = 180; %at apogee true anamoly is 180 deg
h = sqrt(ap*muSun*(1+ecc*cosd(theta)));

%% Calculating Fly By Characteristics

Vsc1 = h/ap;
Ve = sqrt(muSun/Re);
Vinf = abs(Vsc1 - Ve);
phi = 180; %ht ellipse where planet V is faster therefor phi = 180
ecch = 1 + (flyr*(Vinf^2))/mu; %is hyperbolic

%% Finding Turn Angle & phi_2
ta = 2*asind(1/ecch); %turn angle
phi2 = phi - ta;

%% Finding V S/C 2
Vinf2 = [Vinf*cosd(phi2), Vinf*sind(phi2)];
Vsc2 = Vinf2 + [Ve,0];
Vsc1 = [Vsc1, 0];

dV = norm(Vsc1 - Vsc2); %<3: This is under the earth max so its ok
disp(['Gained Delta V (km/s): ', num2str(dV)])

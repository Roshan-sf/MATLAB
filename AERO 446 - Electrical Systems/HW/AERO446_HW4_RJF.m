%% Roshan Jaiswal-Ferri
%Section - 01
%Aero 446 HW3: 4/25/25

%% Workspace Prep

%warning off
format long     %Allows for more accurate decimals
close all;      %Clears all
clear all;      %Clears Workspace
clc;            %Clears Command Window

%% Problem 1

Rearth = 6378; %km
mu = 398600;
p = 6*3600;
beta = 0;
pout = 500;
dod = 0.5;
n = 1;

syms Rorbit
eqn = p == 2*pi*sqrt((Rorbit^3)/mu);
Rorbit = solve(eqn,Rorbit);
Rorbit = double(Rorbit(1));

fe = (1/pi) * asin( sqrt((Rearth/Rorbit)^2 - sin(beta)^2) / cos(beta) );

eclipseTime = 6*3600*fe; %ecplipse time in seconds

energy = (pout*eclipseTime)/(dod*n); %in watt seconds
energy = energy/3600; %convert to watt hours

disp(['Ecplipse Time (s): ', num2str(eclipseTime)])
disp(['Battery size (w-hr): ', num2str(energy)])
disp(' ')

%% Problem 2

Pin = (pout*eclipseTime)/(p-eclipseTime);
area = (Pin)/245;

disp(['Solar Panel Area 1 (m^2): ', num2str(area)])
disp(' ')

%% Problem 3
    
Pin2 = (600*eclipseTime)/(p-eclipseTime);

area2 = Pin2/245;

disp(['Solar Panel Area 2 (m^2): ', num2str(area2)])
disp(['Difference in Area (m^2): ', num2str(area2-area)])
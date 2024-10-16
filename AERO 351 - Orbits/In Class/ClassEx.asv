%Roshan Jaiswal-Ferri
%Section - 01
%Aero 351 In Class Ex: 10/16/24

%% Workspace Prep

format long     %Allows for more accurate decimals
close all;      %Clears all
clear all;      %Clears Workspace
clc;            %Clears Command Window

%% Part 1: Variables

R = [9031.5, -5316.9, -1647.2]; % km
V = [-2.864, 5.1112, -5.0805]; % km/s
mu = 398600; %in km^3/S^2

%% Part2: Calling Function

[h,a,e,nu,i,RAAN,w,p,t,energy] = rv2coes(R,V,mu);

%% Converting Rad to Deg
nu = rad2deg(nu);
RAAN = rad2deg(RAAN);
i = rad2deg(i);
w = rad2deg(w);

%% Displaying Results
disp(['Specific Angular Momentum: ', num2str(h), ' km^2/s'])
disp(['Semi-Major Axis: ', num2str(a), ' km'])
disp(['Eccentricity: ', num2str(e), ' unitless'])
disp(['True Anomaly: ', num2str(nu), ' deg'])
disp(['Inclination: ', num2str(i), ' deg'])
disp(['RAAN: ', num2str(RAAN), ' deg'])
disp(['Argument of Periapsis: ', num2str(w), ' deg'])
disp(['Period: ', num2str((p/3600)), ' hrs'])
disp(['Time since perigee: ', num2str((t/3600)), ' hrs'])
disp(['Specific Energy: ', num2str(energy), ' J/Kg'])










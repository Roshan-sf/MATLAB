%Roshan Jaiswal-Ferri
%Aero 215 HW3 COEs: 11/2/23

%%

close all;      %Clears all
clear all;      %Clears Workspace
clc;            %Clears Command Window

%% Part 1: Variables

R = [-2315.9, 2168.6, 6314.5]; %𝐾̂ [𝑘𝑚]
V = [-3.0599, 6.0645, -3.2044]; %𝐾̂ [𝑘𝑚/𝑠]
mu = 398600; %in km^3/S^2

%% Part2: Calling Function

[a,e,nu,i,RAAN,w,p] = hw3_orbitalCOEs_Jaiswal_ferriRoshan(R,V,mu);

%% Converting Rad to Deg
nu = rad2deg(nu);
RAAN = rad2deg(RAAN);
i = rad2deg(i);
w = rad2deg(w);

%% Displaying Results
disp(['Semi-Major Axis: ', num2str(a), ' km'])
disp(['Eccentricity: ', num2str(e), ' unitless'])
disp(['True Anomaly: ', num2str(nu), ' deg'])
disp(['Inclination: ', num2str(i), ' deg'])
disp(['RAAN: ', num2str(RAAN), ' deg'])
disp(['Argument of Periapsis: ', num2str(w), ' deg'])






%Roshan Jaiswal-Ferri
%Aero 215 HW4 Orbit Transfer: 11/28/23

%%

close all;      %Clears all
clear all;      %Clears Workspace
clc;            %Clears Command Window

%% Defining Variables

R = [6161.56, 454.07, -2537.72]; %𝐾̂ [𝑘𝑚]
V = [0.376, 7.391, 2.224]; %𝐾̂ [𝑘𝑚/𝑠]
mu = 398600; %in km^3/S^2

%% Part 1: Calculating COEs for satelite

[a,e,nu,i,RAAN,w,p] = hw3_orbitalCOEs_Jaiswal_ferriRoshan(R,V,mu);

% Converting Rad to Deg
nu = rad2deg(nu);
RAAN = rad2deg(RAAN);
i = rad2deg(i);
w = rad2deg(w);
p = p/3600; %seconds to hours

% Displaying Results
disp('Results for Initial Orbit:')
disp(['     Eccentricity: ', num2str(e), ' unitless'])
disp(['     Semi-Major Axis: ', num2str(a), ' km'])
disp(['     True Anomaly: ', num2str(nu), ' deg'])
disp(['     Inclination: ', num2str(i), ' deg'])
disp(['     RAAN: ', num2str(RAAN), ' deg'])
disp(['     Argument of Periapsis: ', num2str(w), ' deg'])
disp(['     Period: ', num2str(p), ' hours'])
disp(' ')
disp(' ')

%% Part 2: Geostationary Orbit

%New Vectors:
R_geo = [42157, 0, 0]; %𝐾̂ [𝑘𝑚]
V_geo = [0, 3.07, 0]; %𝐾̂ [𝑘𝑚/𝑠]

[a,e,nu,i,RAAN,w,p] = hw3_orbitalCOEs_Jaiswal_ferriRoshan(R_geo,V_geo,mu);

% Converting Rad to Deg
nu = rad2deg(nu);
RAAN = rad2deg(RAAN);
i = rad2deg(i);
w = rad2deg(w);
p = p/3600; %seconds to hours

% Displaying Results
disp('Results for Final Geostationary Orbit:')
disp(['     Eccentricity: ', num2str(e), ' unitless'])
disp(['     Semi-Major Axis: ', num2str(a), ' km'])
disp(['     True Anomaly: ', num2str(nu), ' deg'])
disp(['     Inclination: ', num2str(i), ' deg'])
disp(['     RAAN: ', num2str(RAAN), ' deg'])
disp(['     Argument of Periapsis: ', num2str(w), ' deg'])
disp(['     Period: ', num2str(p), ' hours'])
















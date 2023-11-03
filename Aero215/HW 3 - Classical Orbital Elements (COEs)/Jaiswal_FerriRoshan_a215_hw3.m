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

[a,e,nu,i,RAAN,w] = hw3_orbitalCOEs_Jaiswal_ferriRoshan(R,V,mu);

disp(num2str(a))
disp(num2str(e))
disp(num2str(nu))
disp(num2str(i))
disp(num2str(RAAN))
disp(num2str(w))






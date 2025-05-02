%% Roshan Jaiswal-Ferri
%Section - 01 
%Aero 303 Thermocouple Plotting: 3/10/25

%% Workspace Prep

format long     %Allows for more accurate decimals
close all;      %Clears all
clear all;      %Clears Workspace
clc;            %Clears Command Window

%% 

w = [0.871,5.053,0.299,0.051,0.009,0.091,2.02,2.143];
a = [24.88888889,35.88888889,6.597222222,1.662222222,0.14625,1.875,...
109.5555556,70.09722222];

figure('Name','Mass vs Cross Sectional Area')
plot(w,a,'.')
xlabel('Weight (g)')
ylabel('Area (cm^2)')
title('Mass vs Cross Sectional Area')
grid on

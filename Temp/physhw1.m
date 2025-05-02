%% Roshan Jaiswal-Ferri
%Section - 45
%PHYS 142 HW1: 4/7/25

%% Workspace Prep

%warning off
format long     %Allows for more accurate decimals
close all;      %Clears all
clear all;      %Clears Workspace
clc;            %Clears Command Window

%% Problem 47

%Data
mass_g = [100, 150, 200, 250];
time_s = [7.8, 9.8, 10.9, 12.4]; 

mass_kg = mass_g / 1000;

%period
T = time_s / 10;
T_squared = T .^ 2;

%new eq: T^2 = (4*pi^2 / k) * m
p = polyfit(mass_kg, T_squared, 1);
slope = p(1);

k = 4 * pi^2 / slope;

disp(['Estimated spring constant k: ', num2str(k), ' N/m']);

%plot
figure;
plot(mass_kg, T_squared, 'o', 'MarkerFaceColor', 'b'); hold on;
plot(mass_kg, polyval(p, mass_kg), 'r-');
xlabel('Mass (kg)');
ylabel('T^2 (s^2)');
title('T^2 vs. Mass');
legend('Data', 'Linear Fit');
grid on;

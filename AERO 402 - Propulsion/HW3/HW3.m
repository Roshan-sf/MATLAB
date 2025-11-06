%% Roshan Jaiswal-Ferri
%Aero 402 Homework 3: 10/20/25

%% Workspace Prep

format long     %Allows for more accurate decimals
close all;      %Clears all
clear all;      %Clears Workspace
clc;            %Clears Command Window

%% Question 3

% for air
gamma = 1.4; % Specific heat ratio for air
R = 1716; % Gas constant [ft*lbf/(lbm*R)]
gc = 32.174; % Gravitational constant [lbm*ft/(lbf*s^2)]
M = linspace(0.1, 2, 200);  % Mach range from 0.1 to 2

MFP = (1/gc) * sqrt(gamma / R) .* M .* (1 + ((gamma - 1)/2) .* M.^2)...
    .^ (-(gamma + 1)/(2*(gamma - 1)));

figure;
plot(M, MFP);
grid on;
xlabel('Mach Number, M');
ylabel('Mass Flow Parameter (MFP)');
title('Mass Flow Parameter vs Mach Number for Air (\gamma = 1.4)');
xlim([0 2]);
ylim([0 max(MFP)*1.1]);

%% Question 4

gamma = 1.4;
R_slug = 1716; % ft*lbf/(slug*R)
gc = 32.174; %lbm*ft/(lbf*s^2)
M_face = 0.5;
mdot_lbm = 102.56;

T_static = 411.69; % Rankine
P_static = 4.36 * 144; % psf
M = 0.85; % Freestream Mach number

T_total = T_static * (1 + ((gamma - 1)/2) * M^2); %isentropic eqs
P_total = P_static * (1 + ((gamma - 1)/2) * M^2)^(gamma/(gamma - 1));

MFP = sqrt(gamma/R_slug) * ( M_face*(1 + (gamma-1)/2*M_face^2)^(-(gamma+1)/(2*(gamma-1))) );
A_mfp = (mdot_lbm*sqrt(T_total)) / (P_total * gc * MFP_rhs)

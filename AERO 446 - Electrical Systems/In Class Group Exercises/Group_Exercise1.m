%% Roshan Jaiswal-Ferri
%Section - 01
%Aero 446 GE1: 4/22/25

%% Workspace Prep

%warning off
format long     %Allows for more accurate decimals
close all;      %Clears all
clear all;      %Clears Workspace
clc;            %Clears Command Window

%% Thermal Equilibrium

%spacecraft deets
abs = 0.6;
ems = 0.4;
ems2 = 1; % sensor array is bb
A = 0.5*0.5*6; %total area
As = 0.5*0.5; %area of 1 side (for sensor area and wetted area)
Am = A-As; %Area minus sensor area

%environment deets
AU = 1.496e+11; %meters
Tsun = 5800;
rsun = 6.9634*10^8; %meters
Asun = 4*pi*rsun^2;
sb = 5.67*10^-8;

%Sun calcs
Qsun = sb*Asun*Tsun^4;
Qi = Qsun/(4*pi*AU^2);
disp(['Incident Solar Heat Flux (W/m^2): ', num2str(Qi)])

%Equi Calcs
syms T
Pabs = abs*Qi*As;
Pemit1 = sb*ems*Am*T^4;
Pemit2 = sb*ems2*As*T^4;

eqn = Pabs == Pemit1 + Pemit2;
T_sol = solve(eqn, T);

%find actual sol (the positive real one)
T_real = double(T_sol);
T_real = T_real(imag(T_real)==0 & T_real>0);
disp(['Equilibrium Temp (K): ', num2str(T_real)])
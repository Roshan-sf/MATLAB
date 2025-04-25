%% Roshan Jaiswal-Ferri
%Section - 01
%Aero 446 HW2: 4/18/25

%% Workspace Prep

%warning off
format long     %Allows for more accurate decimals
close all;      %Clears all
clear all;      %Clears Workspace
clc;            %Clears Command Window

%% Problem 1

%1) Immediate command!
%2) Multi-Step
%3) Prohibited Commands
%4) Multi-Step
%5) Block Command

%% Problem 2

%                   !!!!table at bottom!!!!

%% Problem 3

%Assume you’re designing two electronic subsystems for a spacecraft at NASA.
% The worst-case predicted temperature ranges are
%• System A: -3°C to +45°C
%• System B: -20°C to +60°C
%What are the Flight Acceptance and Prototype qualification temperature 
% limits for these components?

%System A
%Flight Acceptance:
%   -8 to 50 C
%Prototype:
%   -13 to 55 C

%System A
%Flight Acceptance:
%   -25 to 65 C
%Prototype:
%   -30 to 70 C


% Source: https://explorers.larc.nasa.gov/2019APSMEX/MO/pdf_files/gsfc-std-7000a_final_3-28-18.pdf
%Pg 152
%% Problem 4

T = 327;
r = 6050*1000; %radius in meters
sb = 5.67*10^-8;
E = 1; %black body
A = 4*pi*r^2;

Q = E*sb*A*T^4; %W
disp(['Venus Radiated Heat Flux (Watts): ', num2str(Q)])

%% Problem 5

% It is easier to be to cold because it is much easier to generate heat
% than it is to cool off. Also heating up to much can cause more severe
% physical issues (melting etc).

%% Problem 6

AU = 1.496e+11; %meters
Tsun = 5800;
rsun = 6.96*10^8; %meters
Asun = 4*pi*rsun^2;

Qsun = sb*Asun*Tsun^4;

Qi = Qsun/(4*pi*AU^2);
disp(['Incident Solar Heat Flux (W/m^2): ', num2str(Qi)])

%% Problem 7

%Scenario 1) Ideal Radiator (Emissivity = 1)
E = 1;
abs = 1; %absorptivity
A = 1; %m^2
Ar = 2; %radiator area

Pabs = abs*Qi*A;
T = (Pabs/(sb*Ar*E))^(1/4);
disp(['Scenario 1 Temp (K): ', num2str(T)])

%Scenario 2) Emissivity 0f 0.5 and Absorptivity of 0.2
E = 0.5;
abs = 0.2; %absorptivity
A = 1; %m^2
Ar = 2; %radiator area

Pabs = abs*Qi*A;
T = (Pabs/(sb*Ar*E))^(1/4);
disp(['Scenario 2 Temp (K): ', num2str(T)])

%Scenario 3) Increase equilibrium to 10C with previous E and A

T = 283.15; %10C
Ar = Pabs/(sb*E*T^4);
disp(['New Total** Radiator Area (m^2): ', num2str(Ar)])
disp(['Individual Radiator Area (m^2): ', num2str(Ar/2)])
%% Roshan Jaiswal-Ferri
%Section - 01
%Aero 446 HW7: 6/2/25

%% Workspace Prep

%warning off
format long     %Allows for more accurate decimals
close all;      %Clears all
clear all;      %Clears Workspace
clc;            %Clears Command Window

%% Problem 1

% Constants
f = 10e9; %Frequency in Hz
c = 3e8; %Speed of light in m/s
lambda = c / f; % Wavelength in meters
eta = 0.55;
A_baseline = 1; % antenna area m^2
P_baseline = 100; %Tx power w

G_baseline = 10*log10((4*pi*A_baseline*eta)/lambda^2); %gain in dBi

EIRP_baseline = 10*log10(P_baseline) + G_baseline;

P_values = 10:10:200; %Transmit power values (W)
A_values = zeros(size(P_values)); %Antenna area 

% Calculate required antenna area for each transmit power
for i = 1:length(P_values)
    G_required = EIRP_baseline - 10*log10(P_values(i));
    A_values(i) = (lambda^2/(4*pi*eta))*10^(G_required/10);
end

% Plot results
figure;
plot(P_values, A_values, '-o');
grid on;
xlabel('Transmit Power (W)');
ylabel('Required Antenna Area (m^2)');
title('Antenna Area vs. Transmit Power to Maintain Constant EIRP');

%% Problem 2

Re = 6378;
Rgeo = 36000; %km
c = 3e8; %speed light m/s
fc = 4e9; %4Ghz
lambda = c/fc;
b = 50e6; %50 mhz
Ptx = 75; %w
OBO = 3; %dB
LineLoss = 2; %dB
Gr = 21; %dB
F = 3; %dB
phi = 60; %deg
Cn = 10; %dB

z = Re*cosd(phi);
y = Re*sind(phi);
x = Rgeo + Re - z;
FOR = 2*atand(y/x);

beamwidth = deg2rad(FOR);

d = (1.22*lambda)/beamwidth;

disp(['Diameter of Recieving Antenna: ', num2str(d)])

% 10^(dB/10)

%Cn = EIRP + (Gr/Ts) - Ls + 228.6 - (10*log10(B));
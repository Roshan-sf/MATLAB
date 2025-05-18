%% Roshan Jaiswal-Ferri
%Section - 01
%Aero 356: 5/15/25

%% Workspace Prep

%warning off
format long     %Allows for more accurate decimals
close all;      %Clears all
clear all;      %Clears Workspace
clc;            %Clears Command 

%% Question 1

%Rad Weighting
Wp = 2;
We = 1; %weighting for electron and proton

%Body Parts
Ws = 0.01; %skin
Wst = 0.12; %stomach
Wl = 0.12; %lungs
Wb = 0.01; %brain
Bp = [Ws,Wst,Wl,Wb];

Dp = 1.7; %dose protons in grays
De = 1.6; %does electrons in grays

Ht = (Dp*Wp) + (De*We);
E = 0.25*sum(Ht.*Bp); %in sv or Gray they are the same lol
disp(['Absorbed Radiation: ', num2str(E),' Sv'])

%% Question 2

attenuation = 10; %both in mm
x = 2.3;
Tao = exp((-x/attenuation));

abs = -log10(Tao);
disp(['Absorptance: ', num2str(abs)])

%% Question 3

LET = 3e6 * 1e3; % MeV*cm^2*mg^-1 to eV*cm^2*g^-1
Ei = 4.43; %eV
rho = 6.2; %g/cm^3
q = 1.6e-19; %coulumbs

Qd = (LET*rho*q)/Ei;
disp(['Charge per length: ', num2str(Qd), ' C/cm'])

%% Question 4

mu = 0.8; %cm^-1
x = 0.8; %cm

T = 100*exp(-mu*x);
disp(['Percent of passing photons: ', num2str(T)])

%% Question 5

Vp = 6e3; %m/s
B = 0.3e-6; %T
m = 1.6726e-27; % kilograms
q = 1.6e-19; %coulumbs

rl = (m*Vp)/(q*B); %meters
disp(['Larmor Radius: ', num2str(rl), ' m'])
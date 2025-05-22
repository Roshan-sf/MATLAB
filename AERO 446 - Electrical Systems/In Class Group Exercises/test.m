%% Roshan Jaiswal-Ferri
%Section - 01
%Aero 446 GE2: 4/22/25

%% Workspace Prep

%warning off
format long     %Allows for more accurate decimals
close all;      %Clears all
clear all;      %Clears Workspace
clc;            %Clears Command Window

%% S/C & Env Data

mu = 398600;
Re = 6378;
altVec = 300:20:2000;

Rsc = Re + altVec;
beta = 0;
P = 2*pi*sqrt(Rsc.^3/mu);

tEclipse = P.*fe;
figure()
plot(alt, tEclipse)


power.Payload = 138;
power.Structure = 30;
power.Thermal = 80;
power.Power = 25;
power.Comms = 40;
power.Computer = 25;

%% Question 1

fe = (1/pi) * asin( sqrt((Re./Rsc).^2 - sin(beta)^2) / cos(beta) );


tEclipse = P.*fe;
figure()
plot(alt, tEclipse)

%% Question 2

tDay = P - tEclipse;
figure()
plot(alt, tDay)

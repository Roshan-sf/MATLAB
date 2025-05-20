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

power.Payload = 138;
power.Structure = 30;
power.Thermal = 80;
power.Power = 25;
power.Comms = 40;
power.Computer = 25;

%% Question 1

alt = 1111;
r = Re + alt;
beta = deg2rad(30); %assume beta angle of 30 degrees (units are radians)
p = (2*pi)*sqrt((r^3)/(mu)); %period in seconds
ph = p/3600; %period in hrs

inside = ((((Re/r)^2)-(sin(beta)^2))^0.5) / (cos(beta));
fe = (1/pi)*asin(inside); %fraction of orbit spent in eclipse
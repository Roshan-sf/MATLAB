%% Roshan Jaiswal-Ferri
%Section - 01
%Aero 356: 6/6/25

%% Workspace Prep

%warning off
format long     %Allows for more accurate decimals
close all;      %Clears all
clear all;      %Clears Workspace
clc;            %Clears Command 

%% Question 1

TEC = 5.48e16;
c = 2.998e+8;
fsig = 3e+9;

dT = (40.31*TEC)/((c*fsig^2));

dr = dT*c

%% Question 2

ne = 4*10^11;
fcrit = 8.979*sqrt(ne);

fcrit = fcrit*1e-6

%% Queation 3

Re = 6378;
r = 926;
mu = 398600;
vsc = sqrt(mu/(r+Re));
vsc = vsc*1000;
ni = 5*10^11;
Ai = pi*3^2;
q = 1.60217663 * 10^-19;
Ii = q*ni*vsc*Ai;
Ii = Ii*1e3
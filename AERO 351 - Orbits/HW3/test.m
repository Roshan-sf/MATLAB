%Roshan Jaiswal-Ferri
%Section - 01 
%Aero 351 Homework 2: 11/13/24

%% Workspace Prep

format long     %Allows for more accurate decimals
close all;      %Clears all
clear all;      %Clears Workspace
clc;            %Clears Command Window

%% PART 1: 
%for orbit 1
rp1 = 8100;
ra1 = 18900;
muearth = 398600;
h1 = sqrt(2*muearth)*sqrt((rp1*ra1)/(rp1+ra1));
ecc1 = (ra1-rp1)/(ra1+rp1);
a1 = (rp1+ra1)/2;
T1 = ((2*pi)/sqrt(muearth))*a1^(3/2);
%for orbit 1
thetaB = 45;
rb = (h1^2/muearth)*(1/(1+ecc1*cosd(thetaB)));
vazb1 = h1/rb;
vrb1 = (muearth/h1)*ecc1*sind(thetaB);
vb1 = sqrt(vazb1^2+vrb1^2);
flightpath1 = atand(vrb1/vazb1);
%period for phase
thetaC = 150;
Ec = 2*atan(sqrt((1-ecc1)/(1+ecc1))*tand(thetaC/2));
Mc = Ec - ecc1*sin(Ec);
tc = (Mc/(2*pi))*T1;
Eb = 2*atan(sqrt((1-ecc1)/(1+ecc1))*tand(thetaB/2));
Mb = Eb - ecc1*sin(Eb);
tb = (Mb/(2*pi))*T1;
tCB = T1-(tc-tb);
T2 = tCB;
%now orbit 2 info
a2 = (T2*sqrt(muearth)/(2*pi))^(2/3);
ra2 = 2*a2-rb;

h2 = sqrt(2*muearth)*sqrt(rb*ra2/(rb+ra2));
vb2 = h2/rb;
flightpath2 = 0;
%finally, deltav
deltavB = sqrt(vb1^2+vb2^2 - 2*vb1*vb2*cosd(flightpath2-flightpath1));
%need 2*deltavB for the answer
fprintf('--------Problem 6.23-------------------------\n')
fprintf('\n The delta-v is = %g km/s ',2*deltavB)
fprintf('\n---------------------------------------------\n')
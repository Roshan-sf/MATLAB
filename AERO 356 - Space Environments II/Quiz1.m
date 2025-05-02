%% Roshan Jaiswal-Ferri & Stefan Rosu
%Section - 01
%Aero 356 Q1: 4/18/25

%% Workspace Prep

%warning off
format long     %Allows for more accurate decimals
close all;      %Clears all
clear all;      %Clears Workspace
clc;            %Clears Command Window

%%

Re = 6378;
alt = 500;
theta = 30;
r = 0.5;
q = 1366;
alb = 0.4;
As = pi*r^2;

rho = (Re+alt)/Re;
G = (2/3) * ((2*rho + rho^(-2)) - (2 + rho^(-2)) * sqrt(rho^2 - 1));

Q = q*alb*As*G*cosd(theta);

disp(num2str(Q))

%%

r = 0.5;
alt = 885;
Esc = 0.36;
Ee = 0.3;
Et = 255; %K
sb = 5.67*10^-8;
SAe = (4*pi*(Re*1000)^2); %convrt to m
SAsc = 4*pi*r^2;

F = (SAsc/SAe)*(1/2)*(1-(1-(Re/(Re+alt))^2)^0.5);

Qir = Esc*Ee*sb*(Et^4)*SAe*F;

disp(num2str(Qir))
    
%%

r = 3.1/2;
SA = 4*pi*r^2;
Top = 282;
Esc = 0.52;
Ts = 2.75;
sb = 5.6704*10^-8;

Qs = Esc*sb*SA*((Top^4)-(Ts^4));

disp(num2str(Qs))

%%

r = 1.6;
alt = 4311;
SAe = (4*pi*(Re*1000)^2); %convrt to m
SAsc = 4*pi*r^2;

F = (1/2)*(1-(1-(Re/(Re+alt))^2)^0.5);

disp(num2str(F))

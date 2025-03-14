%% Roshan Jaiswal-Ferri
%Section - 02
%Aero 331 HW 3: 3/14/25

%% Workspace Prep

warning off
format long     %Allows for more accurate decimals
close all;      %Clears all
clear all;      %Clears Workspace
clc;            %Clears Command Window

%% Problem 2

A1 = 1e-4; %m^2
A2 = 2e-4; %m^2
A3 = 3e-4; %m^2

b = 0.2; %m
h = 0.4; %m
d = 0.2; %m
G = 26e9; %Pa
t = 0.001; %m
a = b*h;
Vz = 10000; %N

Yc = ((-20)+(3*20))/6; %cm
Iyy = (1*(20^2) + 2*(20^2) + 3*(20^2))*2; %m^4

syms q1 q2 q3 q4 q5 q6 q7

eq1 = q2 == q1 - 2083.33;
eq2 = q3 == q1;
eq3 = q5 == q4 + 6250;
eq4 = q6 == q4;
eq5 = q7 == q6 - q1 - 4166;
eq6 = q4 == 4166 + q7 + q3;
eq7 = Vz*(d+b) == (q1*b)*h + (q5*h)*2*b + (q6*b)*h - (q7*h)*b;
eq8 = (q1*b + q2*h + q3*b - q7*h) == (q6*b + q7*h + q4*b + q5*h);

soln = solve([eq1,eq2,eq3,eq4,eq5,eq6,eq7,eq8],q1,q2,q3,q4,q5,q6,q7);

q1 = double(soln.q1);
q2 = double(soln.q2);
q3 = double(soln.q3);
q4 = double(soln.q4);
q5 = double(soln.q5);
q6 = double(soln.q6);
q7 = double(soln.q7);

a1 = (q1*b + q2*h + q3*b - q7*h)*(2*G*a*t)^(-1); %rads
a2 = (q6*b + q7*h + q4*b - q5*h)*(2*G*a*t)^(-1); %rads
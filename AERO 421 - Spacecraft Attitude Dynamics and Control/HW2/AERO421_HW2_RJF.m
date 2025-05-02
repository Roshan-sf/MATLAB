%% Roshan Jaiswal-Ferri & Stefan Rosu
%Section - 01
%Aero 356 Q1: 4/25/25

%% Workspace Prep

%warning off
format long     %Allows for more accurate decimals
close all;      %Clears all
clear all;      %Clears Workspace
clc;            %Clears Command Window

%%

syms eps_x(t) eps_y(t) s
syms I_x I_y I_z nu

% Define constants
a = (I_z - I_y)/I_x * nu;
b = (I_x - I_z)/I_y * nu;

% Define differential equations
eq1 = diff(eps_x, t) + a * eps_y == 0;
eq2 = diff(eps_y, t) + b * eps_x == 0;

% Take Laplace transforms
L_eq1 = laplace(eq1, t, s);
L_eq2 = laplace(eq2, t, s);

% Display transformed equations
disp('Laplace-transformed equations:')
disp(L_eq1)
disp(L_eq2)







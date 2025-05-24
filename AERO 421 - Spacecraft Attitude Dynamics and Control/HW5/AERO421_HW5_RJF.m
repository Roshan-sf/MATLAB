%% Roshan Jaiswal-Ferri
%Section - 01
%Aero 421: 5/23/25

%% Workspace Prep

%warning off
format long     %Allows for more accurate decimals
close all;      %Clears all
clear all;      %Clears Workspace
clc;            %Clears Command 

%% Initial Conditions

J = diag([1200,2000,2800]);

% Initial condition epsilon_0
epsilon_0 = [0.2; -0.5; 0.3];

% Compute eta_0
eta_0 = sqrt(1 - epsilon_0' * epsilon_0);

% Initial angular velocity vector omega_0
omega_0 = [0.1; -0.05; 0.05]; % rad/s

% Case 1: Initial conditions
epsilon_c1 = [0; 0; 0];
eta_c1 = 1;

% Case 2: Initial conditions
epsilon_c2 = [-0.2; 0.4; 0.2];
eta_c2 = sqrt(1 - epsilon_c2' * epsilon_c2);

tspan = 140;

%% Simulink

out = sim("ADCS_Design_RJF.slx");
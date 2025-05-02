%% Roshan Jaiswal-Ferri
%Section - 01
%AERO 302 Homework 4 - 12/10/24

%% Workspace Prep

format long     %Allows for more accurate decimals
close all;      %Clears all
clear all;      %Clears Workspace
clc;            %Clears Command Window

%% PART 1: Plotting 

%assume flow unity: Re = 1

Re = 1;
theta = linspace(0,2*pi,200);
Cpcyl = 1-4*(sin(theta)).^2;
Cpsph = -(3/Re)*cos(theta);

figure
plot(theta,Cpcyl);
hold on
plot(theta,Cpsph)
xlabel('Radians')
ylabel('Cp')
title('Cp vs Theta')
legend('CP Cylinder', 'CP Sphere', Location="best")
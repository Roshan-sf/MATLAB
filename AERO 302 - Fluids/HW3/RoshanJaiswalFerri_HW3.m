%Roshan Jaiswal-Ferri
%Section - 01
%AERO 302 Homework 3 - 11/20/24

%% Workspace Prep

format long     %Allows for more accurate decimals
close all;      %Clears all
clear all;      %Clears Workspace
clc;            %Clears Command Window

%% PART 1: FK2

T = 293;
theta = 5;
g = -9.81;
y = linspace(0,(5/1000),200);
%muAir = sLawAir(293);
mu = 1.002e-3; %pa*s
rho = 1000; %kg/m^3
h = 5/1000;

zeta = -((rho*g*sind(theta))/mu)*(h-y);

plot(zeta, y);
xlabel('hz')
ylabel('Height mm')
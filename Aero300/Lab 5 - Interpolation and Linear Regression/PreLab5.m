%Roshan Jaiswal-Ferri
%Section - 03
%Aero 300 Lab 4 - Iterative Methods to Solve Matrix Equations: 4/26/24

%% Workspace Prep

format long     %Allows for more accurate decimals
close all;      %Clears all
clear all;      %Clears Workspace
clc;            %Clears Command Window

%% PART 1: 

%Define x and y as coordinate points
x = [-5, -2, 4, 5];
y = [4, -1, 2, -5];

p = polyfit(x,y,3); %using Polyfit  for x and y finding degree 3 coeffecients

disp(p) %displaying results, each coefficient is for their respective degree

% Define x and y domain to plot line graph of polyfit
x1 = -5:0.1:5; %Plotting along x axis -5 to 5
y1 = polyval(p,x1); %using polyval to create y values

% Plot figure
figure
plot(x,y,'x','MarkerEdgeColor','r','MarkerSize',12)
hold on
grid on
plot(x1,y1)
xlabel('X Domain')
ylabel('Y Domain')
title('Line of Best Fit via polyfit Command')

%Polyfit() looks like it uses newton's divided difference because they both
%use interpolation methods to find a line with best fit that uses a
%polynomial of a certain degree. QR factorization could relate to the least
%squares method which finds approx silution for a system to fit a curve.


















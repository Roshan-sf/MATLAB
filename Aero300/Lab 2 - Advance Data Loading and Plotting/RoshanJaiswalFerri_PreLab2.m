%Roshan Jaiswal-Ferri
%Section - 03
%Aero 300 Pre Lab 2 - Advance Data Loading and Plotting: 4/9/24

%%

close all;      %Clears all
clear all;      %Clears Workspace
clc;            %Clears Command Window

%% PART 1: plot()

o = -2*pi; %Setting bounds to variables
p = 2*pi;

theta = linspace(o,p,130); %Creating vector with bounds and stepping

y = pi*sin(theta/2); %Example Function

g = theta/2; %Other Example function

figure; %Creating a figure with overlayed functions using plot
plot(theta, y)
hold on
plot(theta, g) %using the plot command to create two overlaying lines on single figure
grid on;
title('Graph of πSin(θ/2) & θ/2')

%% PART 2: contour()

figure;
x = linspace(-2*pi,2*pi);
y = linspace(0,4*pi);
[X,Y] = meshgrid(x,y);
Z = sin(Y)+cos(X);
contour(X,Y,Z) %the contour function creates a topographic map based on x y coordinates with height z

%% PART 3: surf()

figure;
surf(X,Y,Z) %using the same variables used for the contour command, surf creates an actual 3d graph of x y z data

%% Part 4: streamline()

%load wind
%[startX,startY,startZ] = meshgrid(80,20:10:50,0:5:15);

figure
streamline(X,Y,Z,U,V,W,startX,startY,startZ)
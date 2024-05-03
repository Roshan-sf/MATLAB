%Roshan Jaiswal-Ferri
%Section - 03
%Aero 300 Lab 5 - Iterative Methods to Solve Matrix Equations: 4/26/24

%% Workspace Prep

format long     %Allows for more accurate decimals
close all;      %Clears all
clear all;      %Clears Workspace
clc;            %Clears Command Window

%% PART 1: Plot the Wing - Interpolation

af1 = load("Lab5_Data\Lab5_Data\airfoil_1.txt");
y1 = af1(:,1);
z1 = af1(:,2);

af2 = load("Lab5_Data\Lab5_Data\airfoil_2.txt");
y2 = af2(:,1);
z2 = af2(:,2);

p1 = polyfit(y1,z1,5);
p2 = polyfit(y2,z2,5);

x = linspace(af1(1,1),af1(6,1),100);
y = linspace(af2(1,1),af2(6,1),100);

z3 = linspace(af1(1,2),af1(6,2));
z4 = linspace(af2(1,2),af2(6,2));

x1 = polyval(p1,x);
x2 = polyval(p2,y);

X = [ones(100,1),ones(100,1)];

figure('Name','Wing Plot')
surf(X,[x',y'], [x1',x2'])
hold on
%surf(X,Y, [x1,x2])




%linestyle='none'





%% Roshan Jaiswal-Ferri
%Section - 01 
%Aero 355 Midterm: 2/20/25

%% Workspace Prep

format long     %Allows for more accurate decimals
close all;      %Clears all
clear all;      %Clears Workspace
clc;            %Clears Command Window

%% 
t = 0.005;
h = 0.1 - 0.005;

x = 1;
for L = 0:0.005:0.2
    J1(x) = (4*(L*h)^2)/(((L*2)+(h*2))/t);
    J2(x) = (4*((0.185-L)*h)^2)/((((0.185-L)*2)+(h*2))/t);
    J(x) = J1(x) + J2(x);
    x = x + 1;
end

% L = 0.025:0.005:0.175
%L2 = [0.025:0.005:0.175];
L2 = 0:0.005:0.2;

figure('Name','L vs J')
plot(L2,J)
grid on
xlabel('Length of Left Part')
ylabel('Torsional Constant J (m^4)')
title('L vs J')
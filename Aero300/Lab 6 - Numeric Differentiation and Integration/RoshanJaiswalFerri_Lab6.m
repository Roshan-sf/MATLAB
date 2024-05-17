%Roshan Jaiswal-Ferri
%Section - 03
%Aero 300 Lab 6 - Numeric Differentiation and Integration: 5/16/24

%% Workspace Prep

format long     %Allows for more accurate decimals
close all;      %Clears all
clear all;      %Clears Workspace
clc;            %Clears Command Window

%% PART 1: 

disp = load("dispData.mat");

dt = .0202; %step size
t = disp.t(1:end);
t2 = disp.t(1:end-1);
t3 = disp.t(1:end-2);
f = disp.y;

dsdtF2 = ((f(2:end))-f(1:end-1))/(dt); %two point center diff.
dsdtF3 = ((f(3:end))-f(1:end-2))/(2*dt); %three point center diff.
dsdt2F3 = (f(3:end)-(2*f(2:end-1))+f(1:end-2))/(dt^2);

% figure
% plot(t2,dsdtF2)
% hold on
% grid on
% plot(t,f)
% plot(t3,dsdtF3)

% figure
% subplot(3,1,1)
% plot(t,f,'rx')

%% PART 2:

accel = load("accelData.mat");
accelc = accel.acc;
acc = accelc'; %convert to row
t = accel.t;













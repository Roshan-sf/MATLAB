%% Roshan Jaiswal-Ferri
%Section - 01 
%Aero 355 Midterm: 2/20/25

%% Workspace Prep

format long     %Allows for more accurate decimals
close all;      %Clears all
clear all;      %Clears Workspace
clc;            %Clears Command Window

%% Data Manip

Data = readtable("AERO303_SSWT_DataM.xlsx");

t = Data.DTime_seconds_;
M = Data.Mach;

figure
plot(t,M)
grid on
xlabel('Time (Seconds)')
ylabel('Mach #')
title('Mach vs Time SSWT Run 1')

Data2 = readtable("AERO303_SSWT_Data2M.xlsx");

t2 = Data2.DTime_seconds_;
M2 = Data2.Mach;

figure
plot(t2,M2)
grid on
xlabel('Time (Seconds)')
ylabel('Mach #')
title('Mach vs Time SSWT Run 2')
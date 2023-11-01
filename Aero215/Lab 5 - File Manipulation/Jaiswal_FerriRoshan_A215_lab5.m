%Roshan Jaiswal-Ferri
%Aero 215 Lab 5 - File Manipulation: 10/31/23

%% Header

close all;      %Clears all
clear all;      %Clears Workspace
clc;            %Clears Command Window

%% Part 1: Importing Workspace Data

load("turbineRunData.mat")

%% Part 2: Converting From F to C

T01_C = convtemp(T01_F_,'F','C');
T03_C = convtemp(T03_F_,'F','C');
T04_C = convtemp(T04_F_,'F','C');
T05_C = convtemp(T05_F_,'F','C');

%% Part 2: Converting From PSI to N/m^2 (Pa)

P1_Pa = convpres(P1_psig_,'psi','Pa');
P01_Pa = convpres(P01_psig_,'psi','Pa');
P03_Pa = convpres(P03_psig_,'psi','Pa');
P07_Pa = convpres(P07_psig_,'psi','Pa');

%% Part 3: Graphing

figure;
plot(DTime_seconds_, T01_C,'y')
hold on; 
plot(DTime_seconds_, T03_C,'r')
plot(DTime_seconds_, T04_C,'g')
plot(DTime_seconds_, T05_C,'b')
xlabel('Time (S)');
ylabel('Temperature (C)')
title('Temperature vs Time')
grid on;
legend('T01', 'T03', 'T04', 'T05')

figure;
plot(DTime_seconds_, P1_Pa,'y')
hold on; 
plot(DTime_seconds_, P01_Pa,'r')
plot(DTime_seconds_, P03_Pa,'g')
plot(DTime_seconds_, P07_Pa,'b')
xlabel('Time (S)');
ylabel('Pressure (Pa)')
title('Pressure vs Time')
grid on;
legend('P1', 'P01', 'P03', 'P07')













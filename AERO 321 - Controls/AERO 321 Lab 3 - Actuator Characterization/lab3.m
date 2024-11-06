%Roshan Jaiswal-Ferri
%SLO Prop - TurboJet Mk1 Data Interpretation: 6/24

%% Workspace Prep

format long     %Allows for more accurate decimals
close all;      %Clears all
clear all;      %Clears Workspace
clc;            %Clears Command Window

%% PART 1: 

filename = "gabagoo2.txt";
data = load(filename);

%Time correction
%tC = data(1,1) - 1;

%data(:,1) = (data(:,1) - tC);

%Plotting
figure('Name','PT1')
plot(data(:,1),data(:,2))
grid on
hold on
xlabel('Dist (cm)')
ylabel('Motor Speed (pwm)')
title('Dist v Speed')


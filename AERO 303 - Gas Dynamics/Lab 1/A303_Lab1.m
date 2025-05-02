%% Roshan Jaiswal-Ferri
%Section - 01 
%Aero 303 Lab 1: 01/27/25

%% Workspace Prep

format long     %Allows for more accurate decimals
close all;      %Clears all
clear all;      %Clears Workspace
clc;            %Clears Command Window

%%

data = readtable('run_2_data_group_2.csv');
data2 = readtable("run_3_data_group_2.csv");
time_ms = data{:, 1};
time_ms2 = data2{:,1};
voltage = data{:, 2};
lbf = data{:, 3};
lbf2 = data2{:,3};

% Adjust time to start at 0 and convert to seconds
time_s = (time_ms - time_ms(1)) / 1000;
time_s2 = (time_ms2 - time_ms2(1)) / 1000;

% Plot time in seconds vs pound force
figure;
plot(time_s, lbf, 'LineWidth', 1.5);
xlabel('Time (s)');
ylabel('Force (lbf)');
title('Time vs Force');
grid on;

figure;
plot(time_s2, lbf2, 'LineWidth', 1.5);
xlabel('Time (s)');
ylabel('Force (lbf)');
title('Time vs Force');
grid on;
clc; clear all; close all;

%% TVAC 24-hr Test Data
% Temperature Data taken from a heat/cool 24-hour test with Cal Poly's 
% STELLA (Space Thermal Environment Laboratory for Low-Pressure
% Assessment) on 8/27/2025

% Import Data

timestamps = load("TVAC_24hr_Test_Elapsed_Time.mat");
timestamps = table2array(timestamps.TestLog24HrTest82725);           % hh:mm

temps = load("TVAC_24hr_Test_Temps.mat");
temps = table2array(temps.TestLog24HrTest82725);                     % Celsius

% Fix formatting for timestamp strings
timestamps = string(timestamps);
time = duration(timestamps,'InputFormat','hh:mm');         % minutes:seconds

% Extract avg temp variables for plotting
shroud_avg = temps(:,1);
shroud_max = temps(:,2);
shroud_min = temps(:,3);
platen_avg = temps(:,4);
platen_max = temps(:,5);
platen_min = temps(:,6);

% Plot
figure
plot(time, shroud_avg, time, shroud_max, time, shroud_min, ...
     time, platen_avg, time, platen_max, time, platen_min,"Linewidth", 1);

% Plot data traces
h = plot(time, shroud_avg, time, shroud_max, time, shroud_min, ...
         time, platen_avg, time, platen_max, time, platen_min,"Linewidth", 1);

hold on

% Plot target soak temperatures
soak1 = 61*ones(6);
soak2 = -16*ones(6);
soak3 = 33*ones(8);
soak4 = -6*ones(7);
soak5 = -26*ones(9);
soak6 = -1*ones(7);
soak7 = 25*ones(7);

plot(time(7:12), soak1, "r", "LineWidth", 3, "HandleVisibility","off")
plot(time(20:25), soak2, "r", "LineWidth", 3, "HandleVisibility","off")
plot(time(34:41), soak3, "r", "LineWidth", 3, "HandleVisibility","off")
plot(time(49:55), soak4, "r", "LineWidth", 3, "HandleVisibility","off")
plot(time(61:69), soak5, "r", "LineWidth", 3, "HandleVisibility","off")
plot(time(75:81), soak6, "r", "LineWidth", 3, "HandleVisibility","off")
h_soak = plot(time(86:92), soak7, "r", "LineWidth", 3);   % <-- this one visible

hold off

% Legend: first 6 traces + single soak line
legend([h; h_soak], {"Shroud Avg", "Shroud Max", "Shroud Min", ...
                     "Platen Avg", "Platen Max", "Platen Min", ...
                     "Target Soak"}, "Location", "northeast", "FontSize", 14)

xlabel("Time (hh:mm:ss)", "FontSize", 16)
ylabel("Temperature (°C)", "FontSize", 16)
title("24-Hour Heat/Cool TVAC Test Data", "FontSize", 20)
grid on
grid minor
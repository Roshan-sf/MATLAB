%% Roshan Jaiswal-Ferri
%Summer 2025 TVAC SURP - Dr. A

%% Workspace Prep

warning off
format long     %Allows for more accurate decimals
close all;      %Clears all
%clear all;      %Clears Workspace
clc;            %Clears Command Window

%%

data = readtable("TVAC_24hr_Test_08_27_2025 Wed Aug 27 08_36_45.CSV");

figure('Name','Data')
plot(data.TIME,data.TC01)
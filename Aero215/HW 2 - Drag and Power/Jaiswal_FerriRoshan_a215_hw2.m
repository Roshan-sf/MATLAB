%Roshan Jaiswal-Ferri
%Aero 215 HW 2: 10/19/23

%%

close all;      %Clears all
clear all;      %Clears Workspace
clc;            %Clears Command Window

%% Part 1

%Variables of aircraft parameters
aircraft.WS = 8.27; %in lbf/ft^2
aircraft.AR = 5.63; %unitless
aircraft.CD0 = .00275; %unitless
aircraft.e = .78; %unitless
aircraft.CLM = 1.4; %unitless
aircraft.S = 127; %ft^2

v_max = 70; %in knots

[V, D, P] = hw2_dragPower_Jaiswal_FerriRoshan(v_max, aircraft);
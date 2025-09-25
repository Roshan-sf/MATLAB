%Roshan Jaiswal-Ferri
%Section - 01
%Aero 351 

%% Workspace Prep

format long     %Allows for more accurate decimals
close all;      %Clears all
clear all;      %Clears Workspace
clc;            %Clears Command Window

%% Part 1: 

[r,v,DCM, DCMT] = coes2rv(1,0,50,75,80,398,5);


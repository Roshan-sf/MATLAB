%Roshan Jaiswal-Ferri
%Section - 01
%Aero 351 Homework 1 - Datetime Calcs: 9/23/24

%% Workspace Prep

format long     %Allows for more accurate decimals
close all;      %Clears all
clear all;      %Clears Workspace
clc;            %Clears Command Window

%% PART 1: Converting Mean Solar to JD time

day = 23;
month = 9;
year = 2024;
ut = '12:00:00'; %enter time as a string in hh:mm:ss format (24hr time)



JDtime = tojd(day,month,year,ut);
disp(JDtime)
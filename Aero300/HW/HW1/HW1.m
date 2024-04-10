%Roshan Jaiswal-Ferri
%Section - 03
%Aero 300 Lab 1 - Review of the MATLAB Environment: 4/5/24

%%

close all;      %Clears all
clear all;      %Clears Workspace
clc;            %Clears Command Window

%% PART 1: Chapter 0, Section 4 Question 2

% Part a)
x = .00000001;

y = (tan(x)-x)/x^3;

disp(num2str(y))

%Part b)

clearvars;

x = .000001;

y = (exp(x)+cos(x)-sin(x)-2)/(x^3);

disp(num2str(y))









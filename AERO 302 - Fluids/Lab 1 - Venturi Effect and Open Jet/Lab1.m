%Roshan Jaiswal-Ferri
%Section - 03
%Aero 302 Lab 1 - Venturi Effect and Open Jet: 10/1/24

%% Workspace Prep

format long     %Allows for more accurate decimals
close all;      %Clears all
clear all;      %Clears Workspace
clc;            %Clears Command Window

%% PART 1: Notes

%isentropic until fan, disregard column 6 b/c after fan

%use total pres - 3,4,and 5: rho: 1.225 kg/m^3

%test section 1.0472 m^2\

%Units 1 R = 0.8128 M, .25 R 
A = readmatrix('D3.xlsx');
len = length(A(1,:));


R = .8128;
dy = .25*R;
dx = .25*R;

for i = 1:length(A)
    for j = 1:height(A)

    end
end

%% Surface Plot










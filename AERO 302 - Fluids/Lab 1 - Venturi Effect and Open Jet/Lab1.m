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
% A = readmatrix('D3.xlsx');
% len = length(A(1,:));
% 
% 
% R = .8128;
% dy = .25*R;
% dx = .25*R;
% 
% for i = 1:length(A)
%     for j = 1:height(A)
% 
%     end
% end

%Jet entrainment

%% Surface Plot

% Load data
data = readmatrix('Data\D1.xlsx');
data2 = (readmatrix('Data\D2.xlsx'))'; 
data3 = readmatrix('Data\D3.xlsx'); 
data4 = (readmatrix('Data\D4.xlsx'))'; 
data5 = (readmatrix('Data\D5.xlsx'))'; 
data6 = readmatrix('Data\D6.xlsx'); 

%%


% Generate grid (adjust grid spacing as needed)
[X, Y] = meshgrid(1:size(data, 2), 1:size(data, 1));

filename = [data, data2, data3, data4, data5, data6];

% Create contour plot
figure
contour3(X, Y, data3);
%surf(X, Y, filename);
colorbar;
title('Contour Plot for Data 1'); % Adjust title for each plot
xlabel('X-axis');
ylabel('Y-axis');

figure
contour3(X, Y, data4);
%surf(X, Y, filename);
colorbar;
title('Contour Plot for Data 1'); % Adjust title for each plot
xlabel('X-axis');
ylabel('Y-axis');









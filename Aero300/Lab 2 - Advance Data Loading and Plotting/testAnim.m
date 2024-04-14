%Roshan Jaiswal-Ferri
%Section - 03
%Aero 300 Lab 2 - Advance Data Loading and Plotting: 4/12/24

%%

close all;      %Clears all
clear all;      %Clears Workspace
clc;            %Clears Command Window

%% Testing


T = readtable("Data\Data\convergence.dat");

%disp(T)

x = T.Properties;
y = T.Properties.VariableNames;

edit temp.txt

writecell(y,temp.txt);





























%% --------------Plot 6--------------

% C = load("Data\Data\DENSITY_iteration.mat");
% 
% 
% figure
% for k = 1:799
%     h = C.C{1,k};
%     surf(h)
%     drawnow
% end

%figure
%surf(h)

%fanimator(d)

%drawnow and surfaceplot and set function, keep track of iterations use for
%loop?

% 
% q = input('Continue? y/n');
% 
% 
% if q == 8
%     disp('9')
% else
%     disp('ok')
% end











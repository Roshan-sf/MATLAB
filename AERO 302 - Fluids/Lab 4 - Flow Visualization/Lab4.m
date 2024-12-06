%Roshan Jaiswal-Ferri
%Section - 01 
%Aero 351: Space Debris Removal - 11/13/24

%% Workspace Prep

format long     %Allows for more accurate decimals
close all;      %Clears all
clear all;      %Clears Workspace
clc;            %Clears Command Window

%% PART 1: 

hz = 15;
U1 = 2.429*hz; %speed in mm/s
U = U1/1000; %speed in m/s

t = 54-47; %Time for 5 vorticies
vD6 = 5/t; %vorticies / second

t = 11;
vD5 = 5/t;

t = 12;
vD4 = 5/t;

t = 7;
vD3 = 3/7;

t = 11;
vD2 = 3/11;

t = 16;
vD1 = 3/16;

Di = [3.515,1.897,1.05,0.832,0.626,0.308]; % in inches
D = Di.*0.0254; %diam in meters

f = [vD1,vD2,vD3,vD4,vD5,vD6];

for i = 1:6
    S(i) = (f(i)*D(i))/U; %strouhal number
end

%Re Calc:

rho = 1000;
u = 0.0010016; %dyn visc of water at 20C

for i = 1:6
    Re(i) = (rho*U*D(i))/u;
end


%% HW Calcs

Re2 = [5,20,120,1000,10000];
S2 = [0.666,0.697,0.75,0.857,3.541];

%% Plotting

figure;
hold on;
grid on;

% Plot S vs Re
plot(S, Re, 'bo-', 'LineWidth', 2, 'MarkerSize', 8);

% Plot S2 vs Re2
%plot(S2, Re2, 'r*-', 'LineWidth', 2, 'MarkerSize', 8);

% Add labels, title, and legend
ylabel('Reynolds Number (Re)', 'FontSize', 12);
xlabel('Strouhal Number (S)', 'FontSize', 12);
title('Reynolds Number vs Strouhal Number', 'FontSize', 14);
legend('S vs Re', 'S2 vs Re2', 'Location', 'Best');


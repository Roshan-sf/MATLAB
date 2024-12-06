%Roshan Jaiswal-Ferri
%Section - 01
%AERO 302 Homework 3 - 11/20/24

%% Workspace Prep

format long     %Allows for more accurate decimals
close all;      %Clears all
clear all;      %Clears Workspace
clc;            %Clears Command Window

%% PART 1: FK2

T = 293;
theta = 5;
g = -9.81;
y = linspace(0,(5/1000),200);
mu = 1.002e-3; %pa*s
rho = 1000; %kg/m^3
h = 5/1000;

zeta = -((rho*g*sind(theta))/mu)*(h-y);

plot(zeta, y);
xlabel('hz')
ylabel('Height mm')

%% Cd vs Re

Re = [1,5,20,100,10000,1000000,2000000];
Cd = [11.1628,4.0640,2.0403,1.0834,0.3286,0.34426,0.5357];

figure('Name','Cd vs Re (Linear)')
plot(Re, Cd, 'o', 'MarkerSize', 8, 'LineWidth', 1.5);
xlabel('Re (Linear Scale)');
ylabel('Cd');
title('Cd vs. Re (Linear Scale)');
grid on;
% Add labels to points
for i = 1:length(Re)
    text(Re(i), Cd(i), sprintf('(%g, %.2f)', Re(i), Cd(i)), ...
        'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');
end
title('Cd vs Re (Linear)');

figure('Name','Cd vs Re (Log)')
plot(Re, Cd, 'o', 'MarkerSize', 8, 'LineWidth', 1.5);
set(gca, 'XScale', 'log'); 
xlabel('Re (Logarithmic Scale)');
ylabel('Cd');
title('Cd vs. Re (Logarithmic Scale)');
grid on;
% Add labels to points
for i = 1:length(Re)
    text(Re(i), Cd(i), sprintf('(%g, %.2f)', Re(i), Cd(i)), ...
        'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');
end
title('Cd vs Re (Log)');

%% Workspace Prep

clear all;      %Clears Workspace

%% Strouhal #s

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

Re2 = [5,20,120,1000,3000]; %Calculated form # of vorticies per time in animation
S2 = [0.666,0.697,0.75,0.857,1.341];

%% Plotting

figure;
hold on;
grid on;

% Plot S vs Re
plot(Re, S, 'bo-', 'LineWidth', 2, 'MarkerSize', 8);

% Plot S2 vs Re2
plot(Re2, S2, 'r*-', 'LineWidth', 2, 'MarkerSize', 8);

% Add labels, title, and legend
xlabel('Reynolds Number (Re)', 'FontSize', 12);
ylabel('Strouhal Number (S)', 'FontSize', 12);
title('Reynolds Number vs Strouhal Number', 'FontSize', 14);
legend('Lab 4 Data', 'Ansys Sim', 'Location', 'Best');

%%

xy = readmatrix("CdPos");
theta = linspace(0,360,196);

figure
plot(theta,xy(:,2))
set(gca, 'YDir', 'reverse')
xlabel("Degrees")
ylabel('Cd')
title('Cd Vs Theta')
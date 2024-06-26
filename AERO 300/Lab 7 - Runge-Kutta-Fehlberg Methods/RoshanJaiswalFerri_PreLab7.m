%Roshan Jaiswal-Ferri
%Section - 03
%Aero 300 PreLab 7 - Runge-Kutta-Fehlberg Methods: 5/16/24

%% Workspace Prep

format long     %Allows for more accurate decimals
close all;      %Clears all
clear all;      %Clears Workspace
clc;            %Clears Command Window

%% PART 1: Using ODE45

[t, y] = ode45(@fx, [0 10], [1;0]); %Input function, and plot 0-10 initial value y=1 and y'=0

figure
plot(t,y);
grid on
legend('y','yPrime','Location','best')

function dydt = fx(t,y)
    dydt = [y(2);-(y(1)*y(2))-y(1)]; %isolated y double prime, and create function
end


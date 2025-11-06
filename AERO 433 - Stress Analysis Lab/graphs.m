%% Roshan Jaiswal-Ferri
%Aero 433: 10/30/25

%% Workspace Prep
clc; clear; close all;

%% plotting

load("CanData.mat")
time = linspace(0,84,85);

hoopCl = Clara(:,1);
hoopD = Daniella(:,1);
hoopA = Abby(:,1);
hoopR = Roshan(:,1);

longCl = Clara(:,2);
longD = Daniella(:,2);
longA = Abby(:,2);
longR = Roshan(:,2);

figure('Name','Hoop')
grid on; hold on
plot(time(1:numel(hoopCl)),hoopCl)
plot(time(1:numel(hoopD)),hoopD)
plot(time(1:numel(hoopA)),hoopA)
plot(time(1:numel(hoopR)),hoopR)
legend('Clara','Daniella','Abby','Roshan')

figure('Name','Long')
grid on; hold on
plot(time(1:numel(longCl)),longCl)
plot(time(1:numel(longD)),longD)
plot(time(1:numel(longA)),longA)
plot(time(1:numel(longR)),longR)
legend('Clara','Daniella','Abby','Roshan')
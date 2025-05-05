%% Roshan Jaiswal-Ferri
%Section - 01
%Aero 446 HW1: 4/11/25

%% Workspace Prep

warning off
format long     %Allows for more accurate decimals
close all;      %Clears all
clear all;      %Clears Workspace
clc;            %Clears Command Windowv

%%

heat = readtable("hot.csv");
cool = readtable("CylinderCooling17min.csv");

hot = heat.Temperature1__C_Run_1(1:8000);


avg = [0.0402,0.0250,0.0271]; %slopes for each
mb = 29.519;
mw = 29.453;
mp = 27.763;
m = [mb,mp,mw]/1000; %black, polished, white

l = 1.5 * 0.0254; %inch to meter
d = 0.75 * 0.0254;
A = l*d;

cp = 903; %j/kg-k, can use k because same slope as c

for i = 1:3
    I(i) = (m(i)*cp*avg(i))/A;
end

disp(num2str(I))

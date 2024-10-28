%Roshan Jaiswal-Ferri
%Section - 03
%Aero 302 Lab 2 - Venturi Effect and Open Jet: 10/1/24

%% Workspace Prep

format long     %Allows for more accurate decimals
close all;      %Clears all
clear all;      %Clears Workspace
clc;            %Clears Command Window

%% PART 1: Read Data

strip67s = load('RPM_500\S7G3RPM500STRIP6-7.mat');
strip78s = load("RPM_500\S7G2RPM500STRIP7-8.mat");
strip89s = load("RPM_500\S7G1RPM500STRIP8-9.mat");
strip910s = load("RPM_500\S6G3RPM500STRIP9-10.mat");
strip1011s = load("RPM_500\S6G2RPM500STRIP10-11.mat");
strip1112s = load("RPM_500\S6G1RPM500STRIP11-12.mat");
noStrips = load("RPM_500\S6G1RPM500STRIP.mat");

strip67 = strip67s.P;
strip78 = strip78s.P;
strip89 = strip89s.P;
strip910 = strip910s.P;
strip1011 = strip1011s.P;
strip1112 = strip1112s.P;
noStrip = noStrips.P;

%% Math or stuff idk

avg67 = mean(strip67,1);
avg78 = mean(strip78,1);
avg89 = mean(strip89,1);
avg910 = mean(strip910,1);
avg1011 = mean(strip1011);
avg1112 = mean(strip1112);
avgnoStrp = mean(noStrip,1);

Cp67 = coP(avg67);
Cp78 = coP(avg78);
Cp89 = coP(avg89);
Cp910 = coP(avg910);
Cp1011 = coP(avg1011);
Cp1112 = coP(avg1112);
CpNoS = coP(avgnoStrp);

Cp = [Cp67; Cp78; Cp89; Cp910; Cp1011; Cp1112; CpNoS];

%% Plots:

theta1 = linspace(-180, 0, 12);
theta2 = linspace(0, 180, 12);


figure
for i = 1:7
    plot(theta1,Cp(i,16:27))
    hold on
    grid on
end
for i = 1:7
    plot(theta2,Cp(i,4:15))
    hold on
    grid on
end

set(gca, 'Ydir', 'reverse')
xlabel('Degrees')
ylabel('CP')

%% Finding drag (x component)

diam = 15.875*.01; %cm to m
r = diam/2;
angle = 15; %deg
rad = deg2rad(angle);
arcL = r*rad;

u = 1;
for i = 0:15:180
    Fx(1,u) = -cosd(i);
    u = u + 1;
end

u = 1;
for i = 0:15:180
    Fy(1,u) = sind(i);
    u = u + 1;
end

d = LandD(avg67, arcL, Fx, Fy);

%% Function

function [D,L] = LandD(P,S,Fx,Fy)
    u = 1;
    for i = 4:15
        D(1,u) = P(1,i)*S*Fx(1,i);
        u =+ 1;
    end
    for i = 16:27
        D(1,u) = P(1,i)*S*Fx(1,i);
        u =+ 1;
    end
end

function [Cp] = coP(data)
    for i = 4:27
        Pi = data(1,i);
        Pinf = data(1,2);
        qinf = data(1,1)-data(1,2);

        Cp(1,i) = (Pi-Pinf)/(qinf);
    end
end



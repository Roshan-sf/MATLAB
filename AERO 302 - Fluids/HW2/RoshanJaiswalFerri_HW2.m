%Roshan Jaiswal-Ferri
%Section - 01
%AERO 302 Homework 2 - 11/9/24

%% Workspace Prep

format long     %Allows for more accurate decimals
close all;      %Clears all
clear all;      %Clears Workspace
clc;            %Clears Command Window

%% PART 1: CVA1

Ct = linspace(0,6,600);
V0 = linspace(0,6,600);

Vprop = 0.5*V0.*(sqrt(1+Ct)-1);
eff = V0./(Vprop+V0);

figure('Name','Efficiency V Thrust Coeff')
plot(Ct,eff)
hold on
plot(Ct, Vprop./V0)
xlabel('Thrust Coeff Ct')
ylabel('Efficiency \eta')
legend('\eta', 'Ct')

%% PART 2: Isentropic Flow (ISF1)

a = 343;
G = 1.4;
M = linspace(0,3,200);
P0 = 101325;
T0 = 300;
R0 = 1.225;
Cp = 1.005; %of air @300k
R = 287; %j/Kg-K

for i = 1:length(M)
    [T(1,i), P(1,i), R(1,i)] = IsenCalc2(G, M(i), P0, T0, R0);
end

figure('Name','Isentropic Calculations')
plot(M,(P./P0))
hold on
plot(M,(T./T0))
plot(M,(R./R0))
title('Isentropic Calculations of Air')
xlabel('Mach #')
ylabel('Ratio')
legend('P/P0','T/T0','R/R0',Location="best")
grid on

Lab1D = load('S4G1D3RPM400.mat');
T0 = mean(Lab1D.t);
P0 = mean(Lab1D.P(:,1));
P1 = mean(Lab1D.P(:,2));
P2 = mean(Lab1D.P(:,3));
P3 = mean(Lab1D.P(:,4));
P4 = mean(Lab1D.P(:,5));
P5 = mean(Lab1D.P(:,6));

V1 = sqrt(abs((2*(P0-P1))/R0));
V2 = sqrt(abs((2*(P0-P2))/R0));
V3 = sqrt(abs((2*(P0-P3))/R0));
V4 = sqrt(abs((2*(P0-P4))/R0));
V5 = sqrt(abs((2*(P0-P5))/R0));

M2(1,1) = V1/a;
M2(1,2) = V2/a;
M2(1,3) = V3/a;
M2(1,4) = V4/a;
M2(1,5) = V5/a;

for i = 1:length(M2)
    [T1(1,i), P1(1,i), R1(1,i)] = IsenCalc2(G, M2(i), P0, T0, R0);
end

M2S = sort(M2);
figure('Name','Ratio v Location')
plot(M2,(P1/P0),'*')
hold on
plot(M2,(T1/T0),'*')
plot(M2,(R1/R0),'*')
title('Wind Tunnel Location v Ratio')
xlabel('Wind Tunnel Location')
ylabel('Ratio')
legend('P/P0','T/T0','R/R0',Location="best")
grid on
xticks(M2S)  %Set ticks at each Mach number location
xticklabels({'L1', 'L5', 'L4', 'L2', 'L3'})

dS = Cp*log(T1(1,3)/T1(1,2))-R*log(P1(1,3)/P1(1,2));
dS2 = Cp*log(T1(1,5)/T1(1,4))-R*log(P1(1,5)/P1(1,4));

%% PART 3: ISF2

[T, pt, rho] = statm(10000);

a = 299.5; %Speed of sound at 10 km in m/s

v = linspace(1, 284, 400); %Speed from 1 m/s to 284 m/s
M = v/a;

q = 0.5*rho.*v.^2;
p_static = pt-q;
T_total = T*(1+(G-1)/2.*M.^2); %Total temperature (isentropic relation)

figure('Name','Bernoulli Calcutions');
subplot(3,1,1);
plot(M, p_static);
xlabel('Mach Number');
ylabel('Static Pressure (Pa)');
title('Static Pressure vs Mach Number (Bernoulli)');

subplot(3,1,2);
plot(M, pt * ones(size(M)));
xlabel('Mach Number');
ylabel('Total Pressure (Pa)');
title('Total Pressure vs Mach Number (Bernoulli)');

subplot(3,1,3);
plot(M, T_total);
xlabel('Mach Number');
ylabel('Total Temperature (K)');
title('Total Temperature vs Mach Number (Bernoulli)');

%Isentropic Calcs
p_static_iso = pt./(1+(G-1)/2.*M.^2).^(G/(G-1));
T_static_iso = T./(1+(G-1)/2.*M.^2);
rho_static_iso = rho./(1+(G-1)/2.*M.^2).^(1/(G-1));
T_total_iso = T_static_iso.*(1+(G-1)/2.*M.^2);

figure('Name','Isentropic Calculations');
subplot(3,1,1);
plot(M, p_static_iso);
xlabel('Mach Number');
ylabel('Static Pressure (Pa)');
title('Static Pressure vs Mach Number (Isentropic)');

subplot(3,1,2);
plot(M, pt * ones(size(M)));
xlabel('Mach Number');
ylabel('Total Pressure (Pa)');
title('Total Pressure vs Mach Number (Isentropic)');

subplot(3,1,3);
plot(M, T_total_iso);
xlabel('Mach Number');
ylabel('Total Temperature (K)');
title('Total Temperature vs Mach Number (Isentropic)');


%% Functions

function [T, P, rho] = IsenCalc2(Gamma, M, P0, T0, R0)
    T = ((1+((Gamma-1)/2)*M^2)/T0)^-1;
    P = (((T0/T)^(Gamma/(Gamma-1)))/P0)^-1;
    rho = (((P0/P)^(1/Gamma))/R0)^-1;
end

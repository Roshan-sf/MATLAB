%% Roshan Jaiswal-Ferri
%Section - 01
%Aero 421 HW4: 5/9/25

%% Workspace Prep

%warning off
format long     %Allows for more accurate decimals
close all;      %Clears all
clear all;      %Clears Workspace
clc;            %Clears Command Window

%% Problem 1

JD = 2458981.1666667;

r0 = [5980.8297; -1282.3184; 4125.8019];  % km in ECI
v0 = [1.540887; 7.186813; 0];  % km/s

target_lat = 35.3;       % degrees
target_lon = -120.8;     % degrees
target_alt = 0.2;        % km

target_eci = lla2eci(target_lat, target_lon, target_alt, JD);

%% Functions

function [g, g_norm, T_g] = geo_monte_sim(r_0, v_0, p_hat_0, T_0, jd, errors, n)

end
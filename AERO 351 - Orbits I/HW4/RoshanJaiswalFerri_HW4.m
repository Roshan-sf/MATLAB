%Roshan Jaiswal-Ferri
%Section - 01 
%Aero 351 Homework 4: 12/04/24

%% Workspace Prep

format long     %Allows for more accurate decimals
close all;      %Clears all
clear all;      %Clears Workspace
clc;            %Clears Command Window

%% PART 1: Constant/Global Vars

mu = 398600;
muMars = 42828; %same units as earth (km)
muSun = 1.327e11;
Rearth = 6378; %km
Rmars = 3396; %km
tm = 1;
tol = 1e-8;

%% PART 2: Datetime Calcs

JD2000 = juliandate(2000,1,1,0,0,0);
JD1 = juliandate(2005,8,15,0,0,0);
JD2 = juliandate(2006,3,15,0,0,0);
dt = (JD2 - JD1)*86400; %delta t of transf in seconds

T01 = (JD1-JD2000)/36525;
T02 = (JD2-JD2000)/36525;


%% Earth State Vectors at Departure

PCe = AERO351planetary_elements2(3,T01); %Finds the Planetary Coes of earth
% Mean longitude - longitude of perihelion = Mean Anomaly (L - w_hat = Me)
% w_bar = w + raan
[a,ecc,inc,raan,w,theta,w_hat,L] = pcoes(PCe);


[Re,Ve] = coes2rvd(a,ecc,inc,raan,w,theta,muSun);
Re = Re';
Ve = Ve';

%% Mars State Vectors at Arrival

PCm = AERO351planetary_elements2(4,T02); %Finds the Planetary Coes of mars
% Mean longitude - longitude of perihelion = Mean Anomaly (L - w_hat = Me)
% w_bar = w + raan
[a_m,ecc_m,inc_m,raan_m,w_m,theta_m,w_hat_m,L_m] = pcoes(PCm);


[Rm,Vm] = coes2rvd(a_m,ecc_m,inc_m,raan_m,w_m,theta_m,muSun);
Rm = Rm';
Vm = Vm';
%% Cruise Phase Calcs

[V1,V2] = lambUVBi(Re,Rm,dt,tm,muSun,tol); %starting and ending velocities
%of cruise phase, also V1 & V2 is used to find Vinf 1 & 2 respectively

%% Departure Calcs

Rpark = Rearth + 190;
Vpark = sqrt(mu/Rpark);
%Vinf = abs(norm(V1) - norm(Ve));
Vinf = norm(V1 - Ve);
Vbo = sqrt((Vinf^2)+((2*mu)/Rpark)); %V burn out

dV = abs(Vbo-Vpark); % ~3 km/s

%% Arrival Calcs

Rp = Rmars + 300; %km
Vinf2 = norm(V2 - Vm);
Vbo2 = sqrt((Vinf2^2)+((2*muMars)/Rp)); %V burn out

T = 35*3600; %35hrs in seconds
a = ((muMars*T^2)/(4*pi^2))^(1/3);

%vis-viva eq:
Vp = sqrt(muMars*((2/Rp)-(1/a)));

dV2 = abs(norm(Vbo2) - norm(Vp));
dVt = dV + dV2;

disp(dVt)


%% Planetary Ephemerides from Meeus (1991:202-204) and J2000.0
% Output Vector for Aero351planetary in order:
% planet_coes
% a = semimajor axis (km)
% ecc = eccentricity
% inc = inclination (degrees)
% raan = right ascension of the ascending node (degrees)
% w_hat = longitude of perihelion (degrees)
% L = mean longitude (degrees)

% Inputs:
% planet_id - planet identifier:
% 1 = Mercury
% 2 = Venus
% 3 = Earth
% 4 = Mars
% 5 = Jupiter
% 6 = Saturn
% 7 = Uranus
% 8 = Neptune

function [a,ecc,inc,raan,w,theta,w_hat,L] = pcoes(pcoesVec)
    
    a = pcoesVec(1);
    ecc = pcoesVec(2);
    inc = pcoesVec(3);
    raan = pcoesVec(4);
    w_hat = pcoesVec(5);
    L = pcoesVec(6);

    Me = L - w_hat; %mean anamoly in deg
    Mer = deg2rad(Me); %mean anamoly in rad
    w = w_hat - raan; %arg of peri (deg)

    a2 = sqrt((1-ecc)/(1+ecc)); %not semi major axis
    E = Me2e(Mer,ecc); %eccentric anamoly
    theta = 2*atand(tan(E/2)/a2); %in deg

end
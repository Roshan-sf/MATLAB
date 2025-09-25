%Roshan Jaiswal-Ferri
%Section - 01 
%Aero 351: Space Debris Removal - 11/13/24

%% Workspace Prep

format long     %Allows for more accurate decimals
close all;      %Clears all
clear all;      %Clears Workspace
clc;            %Clears Command Window

%% PART 1: Global Vars

mu = 398600;
Rearth = 6378;

%% PART 2: Reading TLE Data

tleL1 = readTLE("TLE_Data\tle_test.txt");

incL1 = tleL1.inclination;

eccL1 = tleL1.eccentricity;

RAANL1 = tleL1.rightAscension;

%argument of periapsis
ArgPL1 = tleL1.argumentOfPerigee;

%Me given in rev/day:
MeL1 = tleL1.meanAnomaly*((2*pi)/86400); %Me in rad/s

%true anomaly
nuL1 = MetoNu(MeL1,eccL1);

%semi major axis
aL1 = ((mu*MeL1)^(1/3))/MeL1; 

n = tleL1.meanMotion;
n1 = n*((2*pi)/86400);
a = (mu/(n1^2))^(1/3);

%% PART 3: Finding State Vectors (R & V) & Propogation

h = (mu*(aL1*(1-eccL1^2)))^(1/2);
theta = nuL1;
ecc = eccL1;

R = [(((h^2)/mu)*(1/(1+ecc*cosd(theta))))*cosd(theta);...
    (((h^2)/mu)*(1/(1+ecc*cosd(theta))))*sind(theta);0];
V = [(mu/h)*-sind(theta);(mu/h)*(ecc+cosd(theta));0];

[~,Q] = ECI2PERI(ArgPL1,incL1,RAANL1);

R1 = Q'*R;
V1 = Q'*V;

%DEBUG: 
% [~,a,e] = rv2coes(R,V,mu,Rearth);

[~,~,~,~,~,~,~,p] = rv2coes(R,V,mu,Rearth);

PropTime = (JDStart-JDtimeL1)+(5*p); %seconds
timespan = [0, PropTime];
state = [R, V]; 
%INPUTS MUST BE IN THAT ORDER UNTIL OPTIONS
options = odeset('RelTol',1e-8,'AbsTol',1e-8);

[timeNew,stateNew] = ode45(@twobodymotion,timespan,state,options,mu);

RL1 = [stateNew(end,1),stateNew(end,2),stateNew(end,3)];
VL1 = [stateNew(end,4),stateNew(end,5),stateNew(end,6)];

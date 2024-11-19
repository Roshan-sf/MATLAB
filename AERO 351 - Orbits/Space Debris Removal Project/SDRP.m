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

tleL1 = tleread("TLE_Data\LEO1.tle");
tleL2 = tleread("TLE_Data\LEO2.tle");
tleM = tleread("TLE_Data\MEO.tle");
tleG = tleread("TLE_Data\GEO.tle");

incL1 = tleL1.Inclination;
incL2 = tleL2.Inclination;
incM = tleM.Inclination;
incG = tleG.Inclination;

eccL1 = tleL1.Eccentricity;
eccL2 = tleL2.Eccentricity;
eccM = tleM.Eccentricity;
eccG = tleG.Eccentricity;

RAANL1 = tleL1.RightAscensionOfAscendingNode;
RAANL2 = tleL2.RightAscensionOfAscendingNode;
RAANM = tleM.RightAscensionOfAscendingNode;
RAANG = tleG.RightAscensionOfAscendingNode;

%argument of periapsis
ArgPL1 = tleL1.ArgumentOfPeriapsis;
ArgPL2 = tleL1.ArgumentOfPeriapsis;
ArgPM = tleL1.ArgumentOfPeriapsis;
ArgPG = tleL1.ArgumentOfPeriapsis;

%Me given in rev/day:
MeL1 = tleL1.MeanAnomaly*((2*pi)/86400); %Me in rad/s
MeL2 = tleL2.MeanAnomaly*((2*pi)/86400);
MeM = tleM.MeanAnomaly*((2*pi)/86400);
MeG = tleG.MeanAnomaly*((2*pi)/86400);

%true anomaly
nuL1 = MetoNu(MeL1,eccL1);
nuL2 = MetoNu(MeL2,eccL1);
nuM = MetoNu(MeM,eccL1);
nuG = MetoNu(MeG,eccL1);

%semi major axis
aL1 = ((mu*MeL1)^(1/3))/MeL1; 
aL2 = ((mu*MeL2)^(1/3))/MeL2;
aM = ((mu*MeM)^(1/3))/MeM;
aG = ((mu*MeG)^(1/3))/MeG;

%Juilian Date
JDtimeL1 = juliandate(2024,11,13,06,51,24); % CHANGE THESE!!
JDtimeL2 = juliandate(2024,11,13,11,16,53);
JDtimeM = juliandate(2024,11,13,02,18,12);
JDtimeG = juliandate(2024,11,13,11,24,12);
JDStart = juliandate(2024,11,14,0,0,0); %start nov-14-2024 at midnight

%% PART 3: Finding State Vectors (R & V) & Propogation

%INPUT AND OUTPUT IN METERS!!!!!!
[R,V] = keplerian2ijk(aL1*1000,eccL1,incL1,RAANL1,ArgPL1,nuL1);
[R_test,V_test] = keplerian2ijk(aM*1000,eccM,incM,RAANM,ArgPM,nuM);
R = R./1000;
V = V./1000;
R_test = R_test./1000;
V_test = V_test./1000;
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

VL1n = norm(VL1);
Vn = norm(V_test);

dV = 2*VL1n*sind(25);
disp(dV)


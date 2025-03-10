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

%Satellite Selection:
%Leo 1: COURIER 1B
%Leo 2: COURIER 1B ROCKET
%Meo: H2SAT (HEINRICH HERTZ)  
%Geo: GALAXY 14 (G-14)


% Read TLE data using the new readTLE function
tleL1 = readTLE("TLE_Data\LEO1.txt");
tleL2 = readTLE("TLE_Data\LEO2.txt");
tleM = readTLE("TLE_Data\MEO.txt");
tleG = readTLE("TLE_Data\GEO.txt");

% Extract inclination
incL1 = tleL1.inclination;
incL2 = tleL2.inclination;
incM = tleM.inclination;
incG = tleG.inclination;

% Extract eccentricity
eccL1 = tleL1.eccentricity;
eccL2 = tleL2.eccentricity;
eccM = tleM.eccentricity;
eccG = tleG.eccentricity;

% Extract Right Ascension of Ascending Node (RAAN)
RAANL1 = tleL1.rightAscension;
RAANL2 = tleL2.rightAscension;
RAANM = tleM.rightAscension;
RAANG = tleG.rightAscension;

% Extract argument of periapsis
ArgPL1 = tleL1.argumentOfPerigee;
ArgPL2 = tleL2.argumentOfPerigee;
ArgPM = tleM.argumentOfPerigee;
ArgPG = tleG.argumentOfPerigee;

%Me given in rev/day:
MeL1 = tleL1.meanAnomaly*((2*pi)/86400); %Me in rad/s
MeL2 = tleL2.meanAnomaly*((2*pi)/86400);
MeM = tleM.meanAnomaly*((2*pi)/86400);
MeG = tleG.meanAnomaly*((2*pi)/86400);

%true anomaly
nuL1 = MetoNu(MeL1,eccL1);
nuL2 = MetoNu(MeL2,eccL2);
nuM = MetoNu(MeM,eccM);
nuG = MetoNu(MeG,eccG);

%semi major axis
aL1 = Me2a(tleL1.meanMotion,mu);
aL2 = Me2a(tleL2.meanMotion,mu);
aM = Me2a(tleM.meanMotion,mu);
aG = Me2a(tleG.meanMotion,mu);

%Juilian Date
JDtimeL1 = juliandate(2024,11,18,18,01,04); % CHANGE THESE!!
JDtimeL2 = juliandate(2024,11,17,23,04,40);
JDtimeM = juliandate(2024,11,13,02,18,12);
JDtimeG = juliandate(2024,11,13,11,24,12);
JDStart = juliandate(2024,11,22,0,0,0); %start nov-22-2024 at midnight

%% PART 3: Finding State Vectors (R & V) & Propogation

%INPUT AND OUTPUT IN METERS!!!!!!
[R,V] = keplerian2ijk(aL1*1000,eccL1,incL1,RAANL1,ArgPL1,nuL1);
R = R./1000;
V = V./1000;
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


%% PART 4: First Transfer

JDL2 = abs(JDtimeL1-JDtimeL2);
x = 1;
for i = 60:10:4*3600
    dtL = i;
    Ttrnsfr = PropTime+dtL+JDL2;
    tspan = [0,Ttrnsfr];
    [R2,V2] = keplerian2ijk(aL2*1000,eccL2,incL2,RAANL2,ArgPL2,nuL2);
    R2 = R2./1000;
    V2 = V2./1000;
    stateL2 = [R2, V2];
    
    [timeNewL2,stateNewL2] = ode45(@twobodymotion,tspan,stateL2,options,mu);
    RL2 = [stateNewL2(end,1),stateNewL2(end,2),stateNewL2(end,3)];
    VL2 = [stateNewL2(end,4),stateNewL2(end,5),stateNewL2(end,6)];
    
    %parameters for lambert w/ uv function:
    tol = 1e-8;
    tm = 1; %<3: Short way around
    
    [V1,V2] = lambUVBi(RL1,RL2,dtL,tm,mu,tol);
    
    Vf1(x,1:3) = V1 - VL1;
    Vf1n(x) = norm(Vf1(x,1:3));

    Vf2(x,1:3) = VL2 - V2;
    Vf2n(x) = norm(Vf2(x,1:3));

    Vf(x) = Vf1n(x) + Vf2n(x);

    t(x) = i;


    x = x +1;
end
% disp([num2str(Vf)])
% disp([num2str(norm(Vf))])


figure
plot(Vf)
xlabel('Time')
ylabel('\delta V km/s')
title('Required \delta V vs Departure Time')

%%

figure('Name', 'Orbit Trajectory');
plot3(stateNew(:, 1), stateNew(:, 2), stateNew(:, 3), 'b', 'LineWidth', 1.5); % Coasting phase
hold on;
plot3(0, 0, 0, 'g*', 'MarkerSize', 10); % Earth at the origin
plot3(stateNewL2(:, 1), stateNewL2(:, 2), stateNewL2(:, 3), 'r', 'LineWidth', 1.5);
plot3(stateNewL2(end,1),stateNewL2(end,2),stateNewL2(end,3), 'm*', 'MarkerSize', 10); % Earth at the origin
plot3(stateNew(end,1),stateNew(end,2),stateNew(end,3), 'c*', 'MarkerSize', 10); % Earth at the origin
xlabel('X (km)');
ylabel('Y (km)');
zlabel('Z (km)');
grid on;
%legend('Coasting Phase', 'Thrust Orbit', 'Earth');
title('Spacecraft Trajectory under Continuous Thrust');
hold off;

%%

%z = find(Vf == 1.897401728841076);

[z, idx] = min(Vf);

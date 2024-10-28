%Roshan Jaiswal-Ferri
%Section - 01
%Aero 351 Homework 2: 10/16/24

%% Workspace Prep

format long     %Allows for more accurate decimals
close all;      %Clears all
clear all;      %Clears Workspace
clc;            %Clears Command Window

%% PART 1: ODE45 Problem

muEarth = 398600; %km^3/s^2
rVect = [20000 -105000 -19000]; %pos
vVect = [.9 -3.4 -1.5]; %km/s (using km and seconds as units)

timespan = [0 2*60*60]; %5hrs into seconds (units need to match)

state = [rVect vVect]; %creating input variable

%outputs = ode45(@filename w/ eqs of motion,time,state,options,anything else)

%INPUTS MUST BE IN THAT ORDER UNTIL OPTIONS

options = odeset('RelTol',1e-8,'AbsTol',1e-8);%ALWAYS CHANGE THESE

[timeNew,stateNew] = ode45(@twobodymotion,timespan,state,options,muEarth);

endPosVec = [stateNew(end,1), stateNew(end,2), stateNew(end,3)];
origin = [0, 0, 0];

disp('Results for 3.20:')

mag = norm(endPosVec);
disp(['Magnitude of the distance (km): ', num2str(mag)]);


finalVelVec = [stateNew(end,4), stateNew(end,5), stateNew(end,6)];

finalVelMag = norm(finalVelVec);
disp(['Magnitude of the final velocity (km/s): ', num2str(finalVelMag)]);

disp(num2str(endPosVec));
disp(num2str(finalVelVec));
disp(' ')

%% PROBLEM 4.5:

R = [6500, -7500, -2500]; % km
V = [4, 3, -3]; % km/s
mu = 398600; %in km^3/S^2
r = 6378; %radius of earth in km

[h,~,e,nu,i,RAAN,w] = rv2coes(R,V,mu,r);

% Converting Rad to Deg
nu = rad2deg(nu);
RAAN = rad2deg(RAAN);
i = rad2deg(i);
w = rad2deg(w);

% Displaying Results
disp('Results for 4.5:')
disp(['Specific Angular Momentum: ', num2str(h), ' km^2/s'])
disp(['Eccentricity: ', num2str(e), ' unitless'])
disp(['True Anomaly: ', num2str(nu), ' deg'])
disp(['Inclination: ', num2str(i), ' deg'])
disp(['RAAN: ', num2str(RAAN), ' deg'])
disp(['Argument of Periapsis: ', num2str(w), ' deg'])
disp(' ')

%% PROBLEM 4.7:

R = [-6600, -1300, -5200]; % km
V = [-.4, -.5, -.6]; % km/s
mu = 398600; %in km^3/S^2
r = 6378; %radius of earth in km

[~,~,~,~,i] = rv2coes(R,V,mu,r);

% Converting Rad to Deg
i = rad2deg(i);

% Displaying Results
disp('Results for 4.7:')
disp(['Inclination: ', num2str(i), ' deg'])

%% Functions:

function dstate = twobodymotion(time,state,muEarth) %dstate is derivitve of state
    %define vars
    x = state(1);
    y = state(2);
    z = state(3);
    dx = state(4); %vel
    dy = state(5); %vel
    dz = state(6); %vel
    
    %mag of pos vector
    r = norm([x y z]);
    
    %accel: !!eqs of motion!!
    ddx = -muEarth*x/r^3;
    ddy = -muEarth*y/r^3;
    ddz = -muEarth*z/r^3;
    
    dstate = [dx; dy; dz; ddx; ddy; ddz];

end
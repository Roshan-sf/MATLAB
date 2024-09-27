%Roshan Jaiswal-Ferri
%Section - 01
%Aero 351 In Class - Review of ODE45: 9/25/24

%% Workspace Prep

format long     %Allows for more accurate decimals
close all;      %Clears all
clear all;      %Clears Workspace
clc;            %Clears Command Window

%% PART 1: Two Body Equation of Motion

% Vars:
muEarth = 398600; %km^3/s^2
rVect = [3450 -1700 7750];
vVect = [5.4 -5.4 1]; %km/s (using km and seconds as units)

timespan = [0 24*60*60]; %24hrs into seconds (units need to match)

state = [rVect vVect]; %creating input variable

%outputs = ode45(@filename w/ eqs of motion,time,state,options,anything else)

%INPUTS MUST BE IN THAT ORDER UNTIL OPTIONS

options = odeset('RelTol',1e-8,'AbsTol',1e-8);%ALWAYS CHANGE THESE

[timeNew,stateNew] = ode45(@twobodymotion,timespan,state,options,muEarth);

figure('Name','2D Plot')
plot(stateNew(:,1),stateNew(:,2)) %all rows in col 1 then all in col 2
xlabel('x Km')
ylabel('y Km')


figure('Name','3D Plot')
plot3(stateNew(:,1),stateNew(:,2),stateNew(:,3))
xlabel('x Km')
ylabel('y Km')
zlabel('z Km')


%% ODE45 Function
%Follow same order of variables (can have dif names) EXCEPT options

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






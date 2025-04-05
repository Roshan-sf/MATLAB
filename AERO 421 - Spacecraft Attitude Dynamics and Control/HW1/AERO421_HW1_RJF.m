%% Roshan Jaiswal-Ferri
%Section - 01
%Aero 421 HW1: 4/2/25

%% Workspace Prep

%warning off
format long     %Allows for more accurate decimals
close all;      %Clears all
clear all;      %Clears Workspace
clc;            %Clears Command Window

%% Total Mass

scm = 20; %solar cell mass kg
sensorm = 100; %mass of sensor in kg
ms = 500; %kg MehielSat mass

total_mass = (2*scm) + sensorm + ms;

%% Detumble Mode 

%These calculations assume a constant mass distrobution while in the
%detumble mode, disregarding that the components are tucked away in
%different spots

% Center of Mass

com_xd = 0; %com = center of mass
com_yd = 0;
com_zd = 0;

comd = [com_xd, com_yd, com_zd];

% Moment of Inertia
a = 2; %side length in m
Icubed = (1/6)*(total_mass*a^2);

% parrallel axis theorum I + md^2
Ixxd = Icubed;
Iyyd = Icubed;
Izzd = Icubed;

Ixyd = 0;
Ixzd = 0;
Iyzd = 0;

Iyxd = Ixyd;
Izxd = Ixzd;
Izyd = Iyzd;

Id = [Ixxd, Ixyd, Ixzd;...
    Iyxd, Iyyd, Iyzd;...
    Izxd, Izyd, Izzd];

%% Normal Operations

% Center of Mass

com_x = 0; %com = center of mass
com_y = 0;
com_z = (sensorm*(1+0.5))/total_mass;

com = [com_x, com_y, com_z];

% Moment of Inertia
a = 2; %side length in m
Icube = (1/6)*(ms*a^2);

Ipanlx = (1/12)*scm*((3^2)+(0.05^2)); %eq: 1/12*m*(h^2 + w^2)
Ipanly = (1/12)*scm*((0.05^2)+(2^2));
Ipanlz = (1/12)*scm*((3^2)+(2^2));

Isensx = (1/12)*sensorm*((0.25^2)+(1^2));
Isensy = (1/12)*sensorm*((0.25^2)+(1^2));
Isensz = (1/12)*sensorm*((0.25^2)+(0.25^2));

% parrallel axis theorum I + md^2
Ixx = (Icube + ms*com(3)^2) + 2*(Ipanlx + scm*((com(3)^2)+(2.5^2)))...
    + (Isensx + sensorm*(1.5-com(3))^2);
Iyy = (Icube + ms*com(3)^2) + 2*(Ipanly + scm*com(3)^2)...
    + (Isensy + sensorm*(1.5-com(3))^2);
Izz = (Icube) + 2*(Ipanlz + scm*2.5^2) + (Isensz);

Ixy = 0;
Ixz = 0;
Iyz = -2*scm*com(3)*2.5;

Iyx = Ixy;
Izx = Ixz;
Izy = Iyz;

I = [Ixx, Ixy, Ixz;...
    Iyx, Iyy, Iyz;...
    Izx, Izy, Izz];

%% Display Results

disp('Results for Detumbled')
disp(['Total Mass: ', num2str(total_mass)]);
disp(['Center of Mass: ', num2str(comd)]);
disp(' ')
disp(num2str(Id))
disp('---------------------------------------------------------------')

disp('Results for Normal')
disp(['Total Mass: ', num2str(total_mass)]);
disp(['Center of Mass: ', num2str(com)]);
disp(' ')
disp(num2str(I))
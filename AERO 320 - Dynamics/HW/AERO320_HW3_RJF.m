%Roshan Jaiswal-Ferri
%Section - 01
%Aero 320 HW 3 - Spacecraft Problem: 10/15/24

%% Workspace Prep

format long     %Allows for more accurate decimals
close all;      %Clears all
clear all;      %Clears Workspace
clc;            %Clears Command Window

%% PART 1: Finding Rotation Matrix:

rV = [6783; 3391; 1953]; %Position Vector km
vV = [-3.5; 4.39; 4.44]; %Vel Vector km/s

%Converting to F'LVLH

Zlvlh = -(rV/norm(rV));
Ylvlh = -(cross(rV,vV)/norm(cross(rV,vV)));
Xlvlh = cross(Ylvlh,Zlvlh);

%Creating Matrix with new vectors

Clvlh_eci = [Xlvlh, Ylvlh, Zlvlh]';
disp(Clvlh_eci)

%% PART 2: Principle axis of rotation & angle of RM from pt 1

angle = acos((trace(Clvlh_eci)-1)/2); %angle in rads
a1 = (Clvlh_eci(2,3)-Clvlh_eci(3,2))/2/sin(angle);
a2 = (Clvlh_eci(3,1)-Clvlh_eci(1,3))/2/sin(angle);
a3 = (Clvlh_eci(1,2)-Clvlh_eci(2,1))/2/sin(angle);

aV = [a1; a2; a3];

disp(['Angle in Rads: ', num2str(angle)]);
disp(' ')
disp(aV);

%% PART 3: Find the roll, pitch, and yaw angles

roll = atan2(Clvlh_eci(2,3), Clvlh_eci(3,3)); %phi
pitch = -asin(Clvlh_eci(1,3)); %theta
yaw = atan2(Clvlh_eci(1,2), Clvlh_eci(1,1)); %psi

disp(num2str(rad2deg(roll)));
disp(num2str(rad2deg(pitch)));
disp(num2str(rad2deg(yaw)));
disp(' ')

%% PART 4: Find the quaternion associated with Clvlh_eci

eta = sqrt(trace(Clvlh_eci)+1)/2;
eta1 = (Clvlh_eci(2,3)-Clvlh_eci(3,2))/4/(eta);
eta2 = (Clvlh_eci(3,1)-Clvlh_eci(1,3))/4/(eta);
eta3 = (Clvlh_eci(1,2)-Clvlh_eci(2,1))/4/(eta);

epsV = [eta1; eta2; eta3];

fV = [epsV; eta];
disp(fV);


%Roshan Jaiswal-Ferri
%Section - 
%Aero 302 Homework 1 FM1 - 9/26/24

%% Workspace Prep

format long     %Allows for more accurate decimals
close all;      %Clears all
clear all;      %Clears Workspace
clc;            %Clears Command Window

%% PART 1:

%Summary:
%Apollo Capsule Re-entry Velocity-Altitude Map Calculation
%Using numeric integration

%Assumptions:
%Straight down descent angle
%No unexpected natural events

%Vars
g = 9.81; 
BC = 4800; %Ballistic Coefficient (W/(C_D*S)) N/m^2
U0 = 11200; %Initial velocity (m/s) = escape velocity at re-entry
t_final = 500; %Time duration for simulation (seconds)
dt = 1; %Time step (s)

U = zeros(1, t_final/dt);
altitude = zeros(1, t_final/dt);
U(1) = U0;
altitude(1) = 100000; %Initial altitude 100km

%Numerical integration using trapezoid rule
for i = 1:(t_final/dt)
    [~, ~, rho] = stdatm_Jaiswsal_FerriRoshan(altitude(i));
    
    dU_dt = -(1/g) * ((1/BC) * rho * U(i)^2 / 2);
    U(i+1) = U(i) + dt * dU_dt;
    altitude(i+1) = altitude(i) - U(i) * dt; %assuming straight down
   
    if altitude(i+1) < 0 %stop at alt 0
        break;
    end
end

%Plot velocity-altitude map
figure;
plot(U(1:i), altitude(1:i));
xlabel('Velocity (m/s)');
ylabel('Altitude (m)');
title('Velocity vs. Altitude for Apollo Capsule Re-entry');

%Summary:
%Assuming an entry straight down, the velocity decreases with altitude as
%the friction of the atmosphere on apollo slows it down.

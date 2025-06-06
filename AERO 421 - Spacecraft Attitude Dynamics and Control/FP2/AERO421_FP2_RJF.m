%% Roshan Jaiswal-Ferri & Stefan Rosu
%Section - 01
%Aero 421 HW1: 4/2/25
%Note: This was prepared with the help of Dr. Mehiel's Template script

%% Workspace Prep

%warning off
format long     %Allows for more accurate decimals
close all;      %Clears all
clear all;      %Clears Workspace
clc;            %Clears Command Window

%% Part 1 - Mass Properties

J = [812.0396             0             0; %imported from first script
       0      545.3729             0;
       0             0      627.7083];


%% Part 2 - Torque Free Attitudue Simulation

% Spacecraft Orbit Properties (given)
global mu
mu = 398600; % km^3/s^2
h = 53335.2; % km^2/s
e = 0; % none
Omega = 0*pi/180; % radians
inc = 98.43*pi/180; % radians
omega = 0*pi/180; % radians
nu = 0*pi/180; % radians

a = h^2/mu/(1 - e^2);
orbital_period = 2*pi*sqrt(a^3/mu);

% Torque free scenario (Given)
T = [0;0;0];

% Set/Compute initial conditions
% intial orbital position and velocity
[r_ECI_0, v_ECI_0] = coes2rvd(a,e,rad2deg(inc),0,omega,nu,mu);

% Compute inital F_LVLH basis vectors in F_ECI components based on F_LVLH
% definition

rV = r_ECI_0; %Position Vector km
vV = v_ECI_0; %Vel Vector km/s

%Converting to F'LVLH

Zlvlh = -(rV/norm(rV));
Ylvlh = -(cross(rV,vV)/norm(cross(rV,vV)));
Xlvlh = cross(Ylvlh,Zlvlh);

%Creating Matrix with new vectors

Clvlh_eci = [Xlvlh, Ylvlh, Zlvlh]';
disp(num2str(Clvlh_eci))
C_b_ECI_0 = Clvlh_eci;

% Initial Euler angles relating F_body and F_LVLH (given)
phi_0 = 0;
theta_0 = 0;
psi_0 = 0;
E_b_LVLH_0 = [phi_0; theta_0; psi_0];

% Initial Quaternion relating F_body and F_LVLH (given)
q_b_LVLH_0 = [0; 0; 0; 1];

% Compute initial C_LVLH_ECI_0, C_b_LHVL_0, and C_b_ECI_0 rotaiton matrices

% Initial Euler angles relating body to ECI
% E_b_ECI_0 = C2EulerAngles(C_b_ECI_0);
E_b_ECI_0 = rotm2eul(C_b_ECI_0);

% Initial quaternion relating body to E
q_b_ECI_0 = -rotm2quat(C_b_ECI_0);

% Initial body rates of spacecraft (given)
w_b_ECI_0 = [0.001; -0.001; 0.002];

tspan = orbital_period;

out = sim('iOnlyKindOfUnderstandThisRightNow.slx');

%% Plot Results

eul = squeeze(out.E_b_ECI.signals.values);
qua = squeeze(out.q_b_ECI.signals.values);
ang = squeeze(out.w_b_ECI.signals.values);


figure
subplot(3,1,1)
plot(out.tout,ang)
xlabel('time (seconds)')
ylabel('angular velocity (rads/sec)')
title('Angular Velocities')
legend('w_x','w_y','w_z')

subplot(3,1,2)
plot(out.tout,qua)
xlabel('time (seconds)')
ylabel('Quaternion Parameter')
title('Quaternions')
legend('q_1','q_2','q_3','q_4')

subplot(3,1,3)
plot(out.tout,eul)
xlabel('time (seconds)')
ylabel('Angle (deg)')
title('Euler Angles')
legend('\phi','\theta','\psi')

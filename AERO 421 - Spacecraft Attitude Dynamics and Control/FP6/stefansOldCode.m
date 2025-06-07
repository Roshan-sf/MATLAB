clc
close all
clear all
%%
% Mass properties for normal operations phase
zbar = 0.23438;
cm = [0; 0; zbar];
total_mass = 640;
J = [812.0396             0             0
       0      545.3729             0
       0             0      627.7083];

ts = 100; % settling time
zeta = 0.65; % Damping ratio
wn = log(0.02*sqrt(1-zeta^2))/-zeta/ts;
beta = atan(sqrt(1-zeta^2)/zeta);
tr = (pi-beta)/wn/sqrt(1-zeta^2);


syms Mp1

eqn = zeta == sqrt(log(Mp1)^2/(pi^2 + log(Mp1)^2));

Mp = double(solve(eqn, Mp1));

Kp = 2*J*eye(3)*wn^2;
Kd = J*eye(3)*2*zeta*wn;


% eps = [0;0;0];
% eta = 1;
% q_c = [eps; eta];

I_s = diag([0.6, 0.6 ,1.2]);
J = J + I_s;
 
% Spacecraft Orbit Properties (given)
global mu
mu = 398600; % km^3/s^2
h = 53335.2; % km^2/s
e = 0; % none
Omega = 0*pi/180; % radians
inclination = 98.43*pi/180; % radians
omega = 0*pi/180; % radians
nu = 0*pi/180; % radians

a = h^2/mu/(1 - e^2);
orbital_period = 2*pi*sqrt(a^3/mu);

% Torque free scenario (Given)
T = [0;0;0];

% Set/Compute initial conditions
% intial orbital position and velocity
[r_ECI_0, v_ECI_0] = COES2RV_RosuStefan([a,e,rad2deg(inclination),omega,Omega,nu],mu);

% Compute inital F_LVLH basis vectors in F_ECI components based on F_LVLH
% definition

Z_lvlh = -r_ECI_0/norm(r_ECI_0);
Y_lvlh = -cross(r_ECI_0,v_ECI_0)/norm(cross(r_ECI_0,v_ECI_0));
X_lvlh = cross(Y_lvlh,Z_lvlh);

C_LVLH_ECI_0 = [X_lvlh',Y_lvlh',Z_lvlh'];

% Initial Euler angles relating F_body and F_LVLH (given)
phi_0 = 0;
theta_0 = 0;
psi_0 = 0;
E_b_LVLH_0 = [phi_0; theta_0; psi_0];

% Initial Quaternion relating F_body and F_LVLH (given)
q_b_LVLH_0 = [0; 0; 0; 1];

% Compute initial C_LVLH_ECI_0, C_b_LHVL_0, and C_b_ECI_0 rotaiton matrices
%C_LVLH_ECI_0 = [x_LVLH'; y_LVLH'; z_LVLH'];
C_b_LVLH_0 = rotx(rad2deg(phi_0))'*roty(rad2deg(theta_0))'*rotz(rad2deg(psi_0))';
C_b_ECI_0 = C_b_LVLH_0*C_LVLH_ECI_0;

% Initial Euler angles relating body to ECI
E_b_ECI_0 = C2EulerAngles(C_b_ECI_0);

% Initial quaternion relating body to E
q_b_ECI_0 = rotm2quat(C_b_ECI_0);

% Initial body rates of spacecraft (given)
w_b_ECI_0 = [0.001; -0.001; 0.002];

% Set simulation time period
N = 0.2; % Number of Orbits
tspan = orbital_period*N;
tspan = 300;

% Simulate!
out = sim('proj6.slx');

%% Plot Results

subplot(3,1,1)
plot(out.w_b_ECI(:,1), out.w_b_ECI(:,2:4))
title('Angular Velocities')
ylabel('angular velocity (rad/sec)')
legend('\omega_x','\omega_y','\omega_z')
grid on

subplot(3,1,2)
plot(out.tout, out.q_b_ECI(:,2:5))
title('Quaternions')
ylabel('Quaternion Parameter')
legend('q_1','q_2','q_3','q_4')
grid on

subplot(3,1,3)
plot(out.tout, out.E_b_ECI(:,2:4))
title('Euler Angles')
xlabel('time (seconds)')
ylabel('Angle (rad)')
legend('\phi','\theta','\psi')
grid on

sgtitle('Body to ECI Dynamics and Kinematics')
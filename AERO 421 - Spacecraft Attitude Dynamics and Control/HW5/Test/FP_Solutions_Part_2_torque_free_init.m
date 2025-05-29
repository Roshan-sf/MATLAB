%% Final Simulation Project
% Solutions for part 2 - Torque Free Motion
% Aero 421
% Eric Mehiel
% Cal Poly, SLO

clear
close all
clc

%% Part 1 - Mass Properties

% Mass properties for normal operations phase
bus.mass = 500;
sp1.mass = 20;
sp2.mass = 20;
sen.mass = 100;
total_mass = bus.mass + sp1.mass + sp2.mass + sen.mass;

bus.d = 2;
sp1.x = 2;
sp1.y = 3;
sp1.z = 0.05;
sp2.x = 2;
sp2.y = 3;
sp2.z = 0.05;
sen.x = 0.25;
sen.y = 0.25;
sen.z = 1;

bus.r = [0; 0; 0];
sp1.r = [0; -2.5; 0];
sp2.r = [0; 2.5; 0];
sen.r = [0; 0; 1.5];

cm = (bus.mass*bus.r + sp1.mass*sp1.r + sp2.mass*sp2.r + sen.mass*sen.r)/total_mass;
%cm = [0;0;0]
r1 = sp1.r - cm;
r2 = bus.r - cm;
r3 = sp2.r - cm;
r4 = sen.r - cm;

J.sp1 = 1/12*sp1.mass*diag([sp1.y^2 + sp1.z^2; sp1.x^2 + sp1.z^2; sp1.x^2 + sp1.y^2]) - sp1.mass*vcross(r1)*vcross(r1);
J.bus = 1/6*bus.mass*bus.d^2 * eye(3) - bus.mass*vcross(r2)*vcross(r2);
J.sp2 = 1/12*sp2.mass*diag([sp2.y^2 + sp2.z^2; sp2.x^2 + sp2.z^2; sp2.x^2 + sp2.y^2]) - sp2.mass*vcross(r3)*vcross(r3);
J.sensor = 1/12*sen.mass*diag([sen.y^2 + sen.z^2; sen.x^2 + sen.z^2; sen.x^2 + sen.y^2]) - sen.mass*vcross(r4)*vcross(r4);

J.sc = J.sp1 + J.bus + J.sp2 + J.sensor;

J = J.sc;

fprintf('The spacecraft mass for the normal operations mode is:\n')
display(total_mass)
fprintf('The spacecraft center of mass for the normal operations mode is:\n')
display(cm)
fprintf('The Inertia Matrix for the normal operations mode is:\n')
display(J.sc)

%% Part 2 - Torque Free Attitudue Simulation

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
[r_ECI_0, v_ECI_0] = coe2rv([h e Omega inclination omega nu]);

% Compute inital F_LVLH basis vectors in F_ECI components based on F_LVLH
% definition
z_LVLH = -r_ECI_0/norm(r_ECI_0);
y_LVLH = -cross(r_ECI_0, v_ECI_0)/norm(cross(r_ECI_0, v_ECI_0));
x_LVLH = cross(y_LVLH, z_LVLH);

% Initial Euler angles relating F_body and F_LVLH (given)
phi_0 = 0;
theta_0 = 0;
psi_0 = 0;
E_b_LVLH_0 = [phi_0; theta_0; psi_0];

% Initial Quaternion relating F_body and F_LVLH (given)
q_b_LVLH_0 = [0; 0; 0; 1];

% Compute initial C_LVLH_ECI_0, C_b_LHVL_0, and C_b_ECI_0 rotaiton matrices
C_LVLH_ECI_0 = [x_LVLH'; y_LVLH'; z_LVLH'];
C_b_LVLH_0 = Cx(phi_0)*Cy(theta_0)*Cz(psi_0);
C_b_ECI_0 = C_b_LVLH_0*C_LVLH_ECI_0;

% Initial Euler angles relating body to ECI
E_b_ECI_0 = C2EulerAngles(C_b_ECI_0);

% Initial quaternion relating body to E
q_b_ECI_0 = C2quat(C_b_ECI_0);

% Initial body rates of spacecraft (given)
w_b_ECI_0 = [0.001; -0.001; 0.002];

% Set simulation time period
tspan = orbital_period;

% Simulate!
out = sim('test3');

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
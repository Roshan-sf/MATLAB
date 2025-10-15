clc; clear all; close all;

%% Project #1

muearth = 398600; %km^3/s^2
radius_earth = 6378; %km 

%% SDO - Solar Dynamics Observatory
% Take TLE and convert to r and v
inc_tle = 34.1015; %degrees
ecc_tle = 0.0001437;
RAAN_tle = 90.3064; %degrees
w_tle = 125.9101; %degree
Me_tle = 216.3067; %degrees, Mean Anomaly
n_tle = 1.00274146/24/60/60; % rev/sec, mean motion

[r_ECI, v_ECI, ~, ~, T] = tle2rv(n_tle, ecc_tle, Me_tle, inc_tle, RAAN_tle, w_tle);
T_GEO = T; %s
r_initial_GEO = r_ECI;
v_initial_GEO = v_ECI;

%% Visualize satellite orbit
timespan = [0 T_GEO]; %secs, multiplication to get it into seconds
state = [r_initial_GEO; v_initial_GEO];
options = odeset('RelTol',1e-8,'AbsTol',1e-8);
[~,stateGEO_initial] = ode45(@twobodymotion,timespan,state,options,muearth);

% figure
% addEarth();
% hold on
% plot3(stateGEO_initial(:,1),stateGEO_initial(:,2),stateGEO_initial(:,3));
% xlabel("x [km]")
% ylabel("y [km]")
% zlabel("z [km]")
% title("SDO Orbit (Target in GEO)")
% plot3(r_initial_GEO(1,1), r_initial_GEO(2,1), r_initial_GEO(3,1), 'o', 'MarkerSize', 8, 'MarkerFaceColor', 'k')
% plot3(stateGEO_initial(10,1), stateGEO_initial(10,2), stateGEO_initial(10,3), '*', 'MarkerSize', 10, 'Color', 'm')
% legend(" ","Trajectory","Starting Position from TLE","Direction of Travel")
% axis equal 
% grid on

%% Pick starting position of chaser
% Chaser is behind the target


[Q_eci_to_lvlh,~,~] = ECI2RSW(r_initial_GEO,v_initial_GEO);

r_initial_GEO_eci = r_initial_GEO;
v_initial_GEO_eci = v_initial_GEO;

r_initial_chaser_lvlh = [0;-100;0]; % km
v_initial_chaser_lvlh = [0; 0.0050; 0]; % km/s

% LVLH frame spin
h_vec = cross(r_initial_GEO_eci,v_initial_GEO_eci);
h = norm(h_vec);
rmag = norm(r_initial_GEO_eci);
omega_RSW = [0; 0; h/(rmag^2)];      % [rad/s]
rho_eci = Q_eci_to_lvlh'*r_initial_chaser_lvlh;
vrel_eci   = Q_eci_to_lvlh' * (v_initial_chaser_lvlh + cross(omega_RSW, r_initial_chaser_lvlh));% [km/s]

r_chaser_eci = r_initial_GEO + rho_eci;        % absolute chaser position in ECI [km]
v_chaser_eci = v_initial_GEO + vrel_eci;       % absolute chaser velocity in ECI [km/s]


%% Plot initial Target and Chaser locations

% figure
% addEarth();
% hold on
% plot3(stateGEO_initial(:,1),stateGEO_initial(:,2),stateGEO_initial(:,3));
% xlabel("x [km]")
% ylabel("y [km]")
% zlabel("z [km]")
% title("Initial Positions of Target and Chaser")
% plot3(r_initial_GEO(1,1), r_initial_GEO(2,1), r_initial_GEO(3,1), '.', 'MarkerSize', 10, 'MarkerFaceColor', 'k')
% plot3(r_chaser_eci(1,1), r_chaser_eci(2,1), r_chaser_eci(3,1), '*', 'MarkerSize', 10, 'Color', 'm')
% legend(" ","Trajectory","Target","Chaser")
% axis equal 
% grid on


%% Hop 1 code

n = 2*pi/T_GEO; % rad/s, mean motion of target

% --- CHASER RELATIVE (LVLH) in km / km/s ---
% r_initial_chaser_lvlh = r_initial_chaser_lvlh; % km
% v_initial_chaser_lvlh = v_initial_chaser_lvlh; % km/s

% Desired end y in km
% I know where I want to end up
r_end_hop1_chaser_lvlh = [0;40;0]; % km
Delta_y = r_end_hop1_chaser_lvlh(2) - r_initial_chaser_lvlh(2); % +60 km

% One-rev CW relation: Delta_y = -(6*pi/n)*y0dot  =>  y0dot_req = -(n/(6*pi))*Delta_y
y0dot_req   = -(n/(6*pi)) * Delta_y;      % [km/s]
deltaV_hop1 = y0dot_req - v_initial_chaser_lvlh(2);  % [km/s]  (ΔVy to apply at t=0)

% height = 4*y0dot/n; (if needed)

% find hop1 deltaV
t_m1_end = T_GEO;
time = linspace(0,t_m1_end,1000);

% x = -2*y0dot/n*cos(n*t) + 2*y0dot/n;
x_chaser_hop1_lvlh = -2*y0dot_req/n*cos(n*time) + 2*y0dot_req/n;

% propogate total vectors for hop1
r_start_hop1_chaser_lvlh = r_initial_chaser_lvlh; % km
v_start_hop1_chaser_lvlh = v_initial_chaser_lvlh + [0; deltaV_hop1 ;0]; % km/s
r_start_hop1_target_eci = r_initial_GEO_eci;
v_start_hop1_target_eci = v_initial_GEO_eci;

% dydt = mod_linearized_rendezvous(time, state, muearth, n)
% ODE initial state (12x1): [x y z xdot ydot zdot rtx rty rtz vtx vty vtz]
state =  [r_start_hop1_chaser_lvlh;
          v_start_hop1_chaser_lvlh;
          r_start_hop1_target_eci;
          v_start_hop1_target_eci];
[~, new_state] = ode45(@mod_linearized_rendezvous,[0 T_GEO],state,options,muearth,n);

% plot in LVLH (downrange vs altitude)
figure
plot(new_state(:,2),new_state(:,1))
xlabel("downrange -- vbar [km]")
ylabel("altitude -- rbar [km]")
title("Hop 1 from initial chaser location to hold 1 start location")

% find ending velocity of chaser
v_end_hop1_chaser_lvlh = new_state(end,4:6); % km/s
r_end_hop1_target_eci = r_start_hop1_target_eci;
v_end_hop1_target_eci = v_start_hop1_target_eci;

% in order to get plots, need to use dxdotdot, dydotdot, dzdotdot
% put deltaV into those equations and propogate with ode45
% start of CW equations, equations of motion

% i know how far downrange I need to go
% at 100, need to go to 40
% y (downrange distance) = 60
% do the whole project in LVLH first, then go back, convert to eci,
% propogate with ode45, and plot

% %% Hold 1
% current_mission_cost = norm(deltaV_hop1)*2; % km/s
% r_start_hold1_chaser_lvlh = r_end_hop1_chaser_lvlh; % km
% v_start_hold1_chaser_lvlh = new_state(end,4:6) - []; % km/s

%% Football Orbit

Rchaser_LVLH = [0; 40; 0]; %LVLH km
Vchaser_LVLH = [0;	0;	0]; % km/s

Rtarget_ECI = [3.4537e+04; -2.3865e+04; 3.9923e+03]; %ECI km
Vtarget_ECI = [-1.4726; -2.3526; -1.3222]; %km/s

mu = 398600;   % km^3/s^2 (Earth)

%% Chief's mean motion
r = norm(Rtarget_ECI);
v = norm(Vtarget_ECI);
a = 1 / ( 2/r - v^2/mu );      % km, vis-viva
n = sqrt(mu/a^3);              % rad/s (works well for near-circular & short arcs)
P = 2*pi/n;                    % s, one chief orbit


%%

a = 40;
syms dvx
eq1 = a == 2*(dvx/n);
sol = solve(eq1,dvx);
dvx = double(sol);

Vchaser_LVLH_postBurn = [dvx;	0;	0]; % km/s

tspan = [0, P];
state = [Rtarget_ECI; Vtarget_ECI; Rchaser_LVLH; Vchaser_LVLH_postBurn]; 
options = odeset('RelTol',1e-8,'AbsTol',1e-8);

[~,relativeOrbits] = ode45(@relativeMotion,tspan,state,options,mu);

RC = [relativeOrbits(:,7),relativeOrbits(:,8),relativeOrbits(:,9)];
VC = [relativeOrbits(:,10),relativeOrbits(:,11),relativeOrbits(:,12)];

VCnorm = vecnorm(VC, 2, 2);

figure('Name','Relative Distance')
plot(RC(:,2),RC(:,1))
hold on
plot(0, 0, '*')
grid on
xlabel('Y Dist (km)')
ylabel('X Dist (km)')
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

% For all next ECI conversions
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
% figure
% plot(new_state(:,2),new_state(:,1))
% xlabel("downrange -- vbar [km]")
% ylabel("altitude -- rbar [km]")
% title("Hop 1 from initial chaser location to hold 1 start location")

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

new_state_hop1 = new_state;

% ********************** PLOT IN ECI *****************************
% For all next ECI conversions
rho_eci = Q_eci_to_lvlh'*r_start_hop1_chaser_lvlh;
vrel_eci   = Q_eci_to_lvlh' * (v_start_hop1_chaser_lvlh + cross(omega_RSW, r_start_hop1_chaser_lvlh));% [km/s]
r_start_hop1_chaser_eci = r_initial_GEO + rho_eci;        % absolute chaser position in ECI [km]
v_start_hop1_chaser_eci = v_initial_GEO + vrel_eci;       % absolute chaser velocity in ECI [km/s]

timespan = [0 T_GEO]; %secs, multiplication to get it into seconds
state = [r_start_hop1_chaser_eci; v_start_hop1_chaser_eci];
[~,state_hop1_chaser_eci] = ode45(@twobodymotion,timespan,state,options,muearth);

% figure
% addEarth();
% hold on
% plot3(stateGEO_initial(:,1),stateGEO_initial(:,2),stateGEO_initial(:,3),"LineWidth",1.5);
% xlabel("x [km]")
% ylabel("y [km]")
% zlabel("z [km]")
% title("Hop 1 in ECI")
% plot3(state_hop1_chaser_eci(:,1),state_hop1_chaser_eci(:,2),state_hop1_chaser_eci(:,3),"LineWidth",1.5)
% plot3(r_initial_GEO(1,1), r_initial_GEO(2,1), r_initial_GEO(3,1), 'o', 'MarkerSize', 8, 'MarkerFaceColor', 'k')
% plot3(state_hop1_chaser_eci(1,1), state_hop1_chaser_eci(1,2), state_hop1_chaser_eci(1,3), '*', 'MarkerSize', 10, 'Color', 'm')
% plot3(state_hop1_chaser_eci(end,1), state_hop1_chaser_eci(end,2), state_hop1_chaser_eci(end,3), '*', 'MarkerSize', 10, 'Color', 'b')
% legend(" ","Target Trajectory","Chaser Trajectory","Starting Position of Target","Starting position of Chaser","Ending Position of Chaser")
% axis equal 
% grid on

%% Hold 1 - Football
% Football Orbit

Rchaser_LVLH = [0; 40; 0]; %LVLH km
Vchaser_LVLH = [0;	0;	0]; % km/s
r_start_hold1_chaser_lvlh = Rchaser_LVLH;
v_start_hold1_chaser_lvlh = Vchaser_LVLH;

r_start_hold_GEO_eci = r_initial_GEO_eci; %ECI km
v_start_hold_GEO_eci = v_initial_GEO_eci; %km/s

mu = 398600;   % km^3/s^2 (Earth)

% Chief's mean motion
r = norm(r_start_hold_GEO_eci);
v = norm(v_start_hold_GEO_eci);
a = 1 / ( 2/r - v^2/mu );      % km, vis-viva
n = sqrt(mu/a^3);              % rad/s (works well for near-circular & short arcs)
P = 2*pi/n;                    % s, one chief orbit


%

a = 40;
syms dvx
eq1 = a == 2*(dvx/n);
sol = solve(eq1,dvx);
dvx = double(sol);

Vchaser_LVLH_postBurn = [dvx;	0;	0]; % km/s

tspan = [0, P];
state = [r_start_hold_GEO_eci; v_start_hold_GEO_eci; Rchaser_LVLH; Vchaser_LVLH_postBurn]; 
options = odeset('RelTol',1e-8,'AbsTol',1e-8);

[~,relativeOrbits] = ode45(@relativeMotion_Roshan,tspan,state,options,mu);

RC = [relativeOrbits(:,7),relativeOrbits(:,8),relativeOrbits(:,9)];
VC = [relativeOrbits(:,10),relativeOrbits(:,11),relativeOrbits(:,12)];

VCnorm = vecnorm(VC, 2, 2);

% figure('Name','Relative Distance')
% plot(RC(:,2),RC(:,1))
% hold on
% plot(0, 0, '*')
% grid on
% xlabel('Y Dist (km)')
% ylabel('X Dist (km)')

new_state_hold1 = new_state;

% ********************** PLOT IN ECI *****************************
% For all next ECI conversions
rho_eci = Q_eci_to_lvlh'*r_start_hold1_chaser_lvlh;
vrel_eci   = Q_eci_to_lvlh' * (v_start_hold1_chaser_lvlh + cross(omega_RSW, r_start_hold1_chaser_lvlh));% [km/s]
r_start_hold1_chaser_eci = r_initial_GEO + rho_eci;        % absolute chaser position in ECI [km]
v_start_hold1_chaser_eci = v_initial_GEO + vrel_eci;       % absolute chaser velocity in ECI [km/s]

timespan = [0 T_GEO]; %secs, multiplication to get it into seconds
state = [r_start_hold1_chaser_eci; v_start_hold1_chaser_eci];
[~,state_hold1_chaser_eci] = ode45(@twobodymotion,timespan,state,options,muearth);

% figure
% addEarth();
% hold on
% plot3(stateGEO_initial(:,1),stateGEO_initial(:,2),stateGEO_initial(:,3),"LineWidth",1.5);
% xlabel("x [km]")
% ylabel("y [km]")
% zlabel("z [km]")
% title("Hold 1 in ECI - Football")
% plot3(state_hold1_chaser_eci(:,1),state_hold1_chaser_eci(:,2),state_hold1_chaser_eci(:,3),"LineWidth",1.5)
% plot3(r_initial_GEO(1,1), r_initial_GEO(2,1), r_initial_GEO(3,1), 'o', 'MarkerSize', 8, 'MarkerFaceColor', 'k')
% plot3(state_hold1_chaser_eci(1,1), state_hold1_chaser_eci(1,2), state_hold1_chaser_eci(1,3), '*', 'MarkerSize', 14, 'Color', 'm')
% plot3(state_hold1_chaser_eci(end,1), state_hold1_chaser_eci(end,2), state_hold1_chaser_eci(end,3), '.', 'MarkerSize', 14, 'Color', 'b')
% legend(" ","Target Trajectory","Chaser Trajectory","Starting Position of Target","Starting position of Chaser","Ending Position of Chaser")
% axis equal 
% grid on

%% Hop 2 code

% --- CHASER RELATIVE (LVLH) in km / km/s ---
r_start_hop2_chaser_lvlh = RC(end,:)'; % km
v_start_hop2_chaser_lvlh = VC(end,:)'; % km/s
v_start_hop2_chaser_lvlh = [0;v_start_hop2_chaser_lvlh(2);v_start_hop2_chaser_lvlh(3)];

% Desired end y in km
% I know where I want to end up
r_end_hop2_chaser_lvlh = [0;-1;0]; % km
Delta_y = r_end_hop2_chaser_lvlh(2) - r_start_hop2_chaser_lvlh(2); % +60 km

% One-rev CW relation: Delta_y = -(6*pi/n)*y0dot  =>  y0dot_req = -(n/(6*pi))*Delta_y
y0dot_req   = -(n/(6*pi)) * Delta_y;      % [km/s]
deltaV_hop2 = y0dot_req - v_start_hop2_chaser_lvlh(2);  % [km/s]  (ΔVy to apply at t=0)

% height = 4*y0dot/n; (if needed)

% find hop2 deltaV
t_m2_end = T_GEO;
time = linspace(0,t_m2_end,1000);

% x = -2*y0dot/n*cos(n*t) + 2*y0dot/n;
x_chaser_hop2_lvlh = -2*y0dot_req/n*cos(n*time) + 2*y0dot_req/n;

% propogate total vectors for hop1
% r_start_hop2_chaser_lvlh = r_initial_chaser_lvlh; % km
v_start_hop2_chaser_lvlh = v_start_hop2_chaser_lvlh + [0; deltaV_hop2 ;0]; % km/s
r_start_hop2_target_eci = r_initial_GEO_eci;
v_start_hop2_target_eci = v_initial_GEO_eci;

% dydt = mod_linearized_rendezvous(time, state, muearth, n)
% ODE initial state (12x1): [x y z xdot ydot zdot rtx rty rtz vtx vty vtz]
state =  [r_start_hop2_chaser_lvlh;
          v_start_hop2_chaser_lvlh;
          r_start_hop2_target_eci;
          v_start_hop2_target_eci];
[~, new_state] = ode45(@mod_linearized_rendezvous,[0 T_GEO],state,options,muearth,n);

% % plot in LVLH (downrange vs altitude)
% figure
% plot(new_state(:,2),new_state(:,1))
% xlabel("downrange -- vbar [km]")
% ylabel("altitude -- rbar [km]")
% title("Hop 2 from initial chaser location to hold 2 start location (40 km --> 1km)")

% find ending velocity of chaser
v_end_hop2_chaser_lvlh = new_state(end,4:6)'; % km/s
r_end_hop2_target_eci = new_state(7:9)';
v_end_hop2_target_eci = new_state(10:12)';

new_state_hop2 = new_state;

% ********************** PLOT IN ECI *****************************
% For all next ECI conversions
rho_eci = Q_eci_to_lvlh'*r_start_hop2_chaser_lvlh;
vrel_eci   = Q_eci_to_lvlh' * (v_start_hop2_chaser_lvlh + cross(omega_RSW, r_start_hop2_chaser_lvlh));% [km/s]
r_start_hop2_chaser_eci = r_initial_GEO + rho_eci;        % absolute chaser position in ECI [km]
v_start_hop2_chaser_eci = v_initial_GEO + vrel_eci;       % absolute chaser velocity in ECI [km/s]

timespan = [0 T_GEO]; %secs, multiplication to get it into seconds
state = [r_start_hop2_chaser_eci; v_start_hop2_chaser_eci];
[~,state_hop2_chaser_eci] = ode45(@twobodymotion,timespan,state,options,muearth);

% figure
% addEarth();
% hold on
% plot3(stateGEO_initial(:,1),stateGEO_initial(:,2),stateGEO_initial(:,3),"LineWidth",1.5);
% xlabel("x [km]")
% ylabel("y [km]")
% zlabel("z [km]")
% title("Hop 2 in ECI (40 km --> -1 km)")
% plot3(state_hop2_chaser_eci(:,1),state_hop2_chaser_eci(:,2),state_hop2_chaser_eci(:,3),"LineWidth",1.5)
% plot3(r_initial_GEO(1,1), r_initial_GEO(2,1), r_initial_GEO(3,1), 'o', 'MarkerSize', 8, 'MarkerFaceColor', 'k')
% plot3(state_hop2_chaser_eci(1,1), state_hop2_chaser_eci(1,2), state_hop2_chaser_eci(1,3), '*', 'MarkerSize', 14, 'Color', 'm')
% plot3(state_hop2_chaser_eci(end,1), state_hop2_chaser_eci(end,2), state_hop2_chaser_eci(end,3), '*', 'MarkerSize', 14, 'Color', 'b')
% legend(" ","Target Trajectory","Chaser Trajectory","Starting Position of Target","Starting position of Chaser","Ending Position of Chaser")
% axis equal 
% grid on

%% Hold 2 - Oscillatory

r_start_hold2_chaser_lvlh = r_end_hop2_chaser_lvlh; %LVLH km
v_start_hold2_chaser_lvlh = v_end_hop2_chaser_lvlh + [0;-v_end_hop2_chaser_lvlh(2);0];
v_start_hold2_chaser_lvlh = [0;0;0];
deltaV_counterhop = -v_end_hop2_chaser_lvlh(2);

t = T_GEO/4;
z = 1; %km
syms zdot0
eq2 = z == (zdot0/n)*sin(n*t);
sol = solve(eq2,zdot0);
zdot0 = double(sol);

v_start_hold2_chaser_lvlh = v_start_hold2_chaser_lvlh + [0; 0; zdot0];

tspan = [0, T_GEO];
state = [r_initial_GEO_eci; v_initial_GEO_eci; r_start_hold2_chaser_lvlh; v_start_hold2_chaser_lvlh]; 
options = odeset('RelTol',1e-8,'AbsTol',1e-8);

[~,relativeOrbits] = ode45(@relativeMotion_Roshan,tspan,state,options,mu);

RC = [relativeOrbits(:,7),relativeOrbits(:,8),relativeOrbits(:,9)];
VC = [relativeOrbits(:,10),relativeOrbits(:,11),relativeOrbits(:,12)];

% figure('Name','Oscillate')
% plot3(RC(:,2),RC(:,1),RC(:,3))
% hold on
% plot3(0, 0, 0, '*')
% grid on
% xlabel('Y Dist (km)')
% ylabel('X Dist (km)')

new_state_hold2 = new_state;

% ********************** PLOT IN ECI *****************************
% For all next ECI conversions
rho_eci = Q_eci_to_lvlh'*r_start_hold2_chaser_lvlh;
vrel_eci   = Q_eci_to_lvlh' * (v_start_hold2_chaser_lvlh + cross(omega_RSW, r_start_hold2_chaser_lvlh));% [km/s]
r_start_hold2_chaser_eci = r_initial_GEO + rho_eci;        % absolute chaser position in ECI [km]
v_start_hold2_chaser_eci = v_initial_GEO + vrel_eci;       % absolute chaser velocity in ECI [km/s]

timespan = [0 T_GEO]; %secs, multiplication to get it into seconds
state = [r_start_hold2_chaser_eci; v_start_hold2_chaser_eci];
[~,state_hold2_chaser_eci] = ode45(@twobodymotion,timespan,state,options,muearth);

% figure
% addEarth();
% hold on
% plot3(stateGEO_initial(:,1),stateGEO_initial(:,2),stateGEO_initial(:,3),"LineWidth",1.5);
% xlabel("x [km]")
% ylabel("y [km]")
% zlabel("z [km]")
% title("Hold 2 in ECI - Oscillatory")
% plot3(state_hold2_chaser_eci(:,1),state_hold2_chaser_eci(:,2),state_hold2_chaser_eci(:,3),"LineWidth",1.5)
% plot3(r_initial_GEO(1,1), r_initial_GEO(2,1), r_initial_GEO(3,1), 'o', 'MarkerSize', 8, 'MarkerFaceColor', 'k')
% plot3(state_hold2_chaser_eci(1,1), state_hold2_chaser_eci(1,2), state_hold2_chaser_eci(1,3), '*', 'MarkerSize', 14, 'Color', 'm')
% plot3(state_hold2_chaser_eci(end,1), state_hold2_chaser_eci(end,2), state_hold2_chaser_eci(end,3), '.', 'MarkerSize', 14, 'Color', 'g')
% legend(" ","Target Trajectory","Chaser Trajectory","Starting Position of Target","Starting position of Chaser","Ending Position of Chaser")
% axis equal 
% grid on

%% Maneuver 3: Hop
% --- CHASER RELATIVE (LVLH) in km / km/s ---
r_start_hop3_chaser_lvlh = RC(end,:)'; % km
v_start_hop3_chaser_lvlh = VC(end,:)'; % km/s
v_start_hop3_chaser_lvlh = [0;v_start_hop3_chaser_lvlh(2);v_start_hop3_chaser_lvlh(3)];

% Desired end y in km
% I know where I want to end up
r_end_hop3_chaser_lvlh = [0;-300/1000;0]; % km
Delta_y = r_end_hop3_chaser_lvlh(2) - r_start_hop3_chaser_lvlh(2); % +60 km

% One-rev CW relation: Delta_y = -(6*pi/n)*y0dot  =>  y0dot_req = -(n/(6*pi))*Delta_y
y0dot_req   = -(n/(6*pi)) * Delta_y;      % [km/s]
deltaV_hop3 = y0dot_req - v_start_hop3_chaser_lvlh(2);  % [km/s]  (ΔVy to apply at t=0)

% height = 4*y0dot/n; (if needed)

% find hop2 deltaV
t_m3_end = T_GEO;
time = linspace(0,t_m3_end,1000);

% x = -2*y0dot/n*cos(n*t) + 2*y0dot/n;
x_chaser_hop3_lvlh = -2*y0dot_req/n*cos(n*time) + 2*y0dot_req/n;

% propogate total vectors for hop1
% r_start_hop2_chaser_lvlh = r_initial_chaser_lvlh; % km
v_start_hop3_chaser_lvlh = [v_start_hop3_chaser_lvlh(1);v_start_hop3_chaser_lvlh(2);0] + [0; deltaV_hop3 ;0]; % km/s
r_start_hop3_target_eci = r_initial_GEO_eci;
v_start_hop3_target_eci = v_initial_GEO_eci;

% dydt = mod_linearized_rendezvous(time, state, muearth, n)
% ODE initial state (12x1): [x y z xdot ydot zdot rtx rty rtz vtx vty vtz]
state =  [r_start_hop3_chaser_lvlh;
          v_start_hop3_chaser_lvlh;
          r_start_hop3_target_eci;
          v_start_hop3_target_eci];
[~, new_state] = ode45(@mod_linearized_rendezvous,[0 T_GEO],state,options,muearth,n);

% % plot in LVLH (downrange vs altitude)
% figure
% plot(new_state(:,2)*1000,new_state(:,1)*1000)
% xlabel("downrange -- vbar [m]")
% ylabel("altitude -- rbar [m]")
% title("Hop 3 from -1 km vbar chaser location to hold 3 start location (-1 km --> -300 m)")

% find ending velocity of chaser
v_end_hop3_chaser_lvlh = new_state(end,4:6)'; % km/s
r_end_hop3_target_eci = new_state(7:9)';
v_end_hop3_target_eci = new_state(10:12)';

new_state_hop3 = new_state;

% ********************** PLOT IN ECI *****************************
% For all next ECI conversions
rho_eci = Q_eci_to_lvlh'*r_start_hop3_chaser_lvlh;
vrel_eci   = Q_eci_to_lvlh' * (v_start_hop3_chaser_lvlh + cross(omega_RSW, r_start_hop3_chaser_lvlh));% [km/s]
r_start_hop3_chaser_eci = r_initial_GEO + rho_eci;        % absolute chaser position in ECI [km]
v_start_hop3_chaser_eci = v_initial_GEO + vrel_eci;       % absolute chaser velocity in ECI [km/s]

timespan = [0 T_GEO]; %secs, multiplication to get it into seconds
state = [r_start_hop3_chaser_eci; v_start_hop3_chaser_eci];
[~,state_hop3_chaser_eci] = ode45(@twobodymotion,timespan,state,options,muearth);

% figure
% addEarth();
% hold on
% plot3(stateGEO_initial(:,1),stateGEO_initial(:,2),stateGEO_initial(:,3),"LineWidth",1.5);
% xlabel("x [km]")
% ylabel("y [km]")
% zlabel("z [km]")
% title("Hop 3 in ECI (-1 km --> -300 m)")
% plot3(state_hop3_chaser_eci(:,1),state_hop3_chaser_eci(:,2),state_hop3_chaser_eci(:,3),"LineWidth",1.5)
% plot3(r_initial_GEO(1,1), r_initial_GEO(2,1), r_initial_GEO(3,1), 'o', 'MarkerSize', 8, 'MarkerFaceColor', 'k')
% plot3(state_hop3_chaser_eci(1,1), state_hop3_chaser_eci(1,2), state_hop3_chaser_eci(1,3), '*', 'MarkerSize', 14, 'Color', 'm')
% plot3(state_hop3_chaser_eci(end,1), state_hop3_chaser_eci(end,2), state_hop3_chaser_eci(end,3), '*', 'MarkerSize', 14, 'Color', 'b')
% legend(" ","Target Trajectory","Chaser Trajectory","Starting Position of Target","Starting position of Chaser","Ending Position of Chaser")
% axis equal 
% grid on

%% Hold 3: V-bar Station Keeping

r_start_hold3_chaser_lvlh = r_end_hop3_chaser_lvlh; %LVLH km
v_start_hold3_chaser_lvlh = [0;	0; 0]; % km/s

tspan = [0, T_GEO];
state = [r_initial_GEO_eci; v_initial_GEO_eci; r_start_hold3_chaser_lvlh; v_start_hold3_chaser_lvlh]; 
options = odeset('RelTol',1e-8,'AbsTol',1e-8);

[~,relativeOrbits] = ode45(@relativeMotion_Roshan,tspan,state,options,mu);

RC = [relativeOrbits(:,7),relativeOrbits(:,8),relativeOrbits(:,9)];
VC = [relativeOrbits(:,10),relativeOrbits(:,11),relativeOrbits(:,12)];

VCnorm = vecnorm(VC, 2, 2);

% figure('Name','Relative Distance')
% plot(RC(:,2)*1000,RC(:,1)*1000)
% hold on
% plot(0, 0, '*')
% grid on
% xlabel('Y Dist (m)')
% ylabel('X Dist (m)')

new_state_hold3 = new_state;

% ********************** PLOT IN ECI *****************************
% For all next ECI conversions
rho_eci = Q_eci_to_lvlh'*r_start_hold3_chaser_lvlh;
vrel_eci   = Q_eci_to_lvlh' * (v_start_hold3_chaser_lvlh + cross(omega_RSW, r_start_hold3_chaser_lvlh));% [km/s]
r_start_hold3_chaser_eci = r_initial_GEO + rho_eci;        % absolute chaser position in ECI [km]
v_start_hold3_chaser_eci = v_initial_GEO + vrel_eci;       % absolute chaser velocity in ECI [km/s]

timespan = [0 T_GEO]; %secs, multiplication to get it into seconds
state = [r_start_hold3_chaser_eci; v_start_hold3_chaser_eci];
[~,state_hold3_chaser_eci] = ode45(@twobodymotion,timespan,state,options,muearth);

% figure
% addEarth();
% hold on
% plot3(stateGEO_initial(:,1),stateGEO_initial(:,2),stateGEO_initial(:,3),"--","LineWidth",1.5);
% xlabel("x [km]")
% ylabel("y [km]")
% zlabel("z [km]")
% title("Hold 3 in ECI - V-Bar Station Keeping at -300 m relative to Targt")
% plot3(state_hold3_chaser_eci(:,1),state_hold3_chaser_eci(:,2),state_hold3_chaser_eci(:,3),"LineWidth",1)
% plot3(r_initial_GEO(1,1), r_initial_GEO(2,1), r_initial_GEO(3,1), 'o', 'MarkerSize', 8, 'MarkerFaceColor', 'k')
% plot3(state_hold3_chaser_eci(1,1), state_hold3_chaser_eci(1,2), state_hold3_chaser_eci(1,3), '*', 'MarkerSize', 14, 'Color', 'm')
% plot3(state_hold3_chaser_eci(end,1), state_hold3_chaser_eci(end,2), state_hold3_chaser_eci(end,3), '.', 'MarkerSize', 14, 'Color', 'g')
% legend(" ","Target Trajectory","Chaser Trajectory","Starting Position of Target","Starting position of Chaser","Ending Position of Chaser")
% axis equal 
% grid on

%% Maneuver 4: Hop

% --- CHASER RELATIVE (LVLH) in km / km/s ---
r_start_hop4_chaser_lvlh = RC(end,:)'; % km
v_start_hop4_chaser_lvlh = VC(end,:)'; % km/s
v_start_hop4_chaser_lvlh = [0;v_start_hop4_chaser_lvlh(2);v_start_hop4_chaser_lvlh(3)];

% Desired end y in km
% I know where I want to end up
r_end_hop4_chaser_lvlh = [0;-20/1000;0]; % km
Delta_y = r_end_hop4_chaser_lvlh(2) - r_start_hop4_chaser_lvlh(2); % +60 km

% One-rev CW relation: Delta_y = -(6*pi/n)*y0dot  =>  y0dot_req = -(n/(6*pi))*Delta_y
y0dot_req   = -(n/(6*pi)) * Delta_y;      % [km/s]
deltaV_hop4 = y0dot_req - v_start_hop4_chaser_lvlh(2);  % [km/s]  (ΔVy to apply at t=0)

% height = 4*y0dot/n; (if needed)

% find hop2 deltaV
t_m4_end = T_GEO;
time = linspace(0,t_m4_end,1000);

% x = -2*y0dot/n*cos(n*t) + 2*y0dot/n;
x_chaser_hop4_lvlh = -2*y0dot_req/n*cos(n*time) + 2*y0dot_req/n;

% propogate total vectors for hop1
% r_start_hop2_chaser_lvlh = r_initial_chaser_lvlh; % km
v_start_hop4_chaser_lvlh = v_start_hop4_chaser_lvlh + [0; deltaV_hop4 ;0]; % km/s
r_start_hop4_target_eci = r_initial_GEO_eci;
v_start_hop4_target_eci = v_initial_GEO_eci;

% dydt = mod_linearized_rendezvous(time, state, muearth, n)
% ODE initial state (12x1): [x y z xdot ydot zdot rtx rty rtz vtx vty vtz]
state =  [r_start_hop4_chaser_lvlh;
          v_start_hop4_chaser_lvlh;
          r_start_hop4_target_eci;
          v_start_hop4_target_eci];
[~, new_state] = ode45(@mod_linearized_rendezvous,[0 T_GEO],state,options,muearth,n);

% % plot in LVLH (downrange vs altitude)
% figure
% plot(new_state(:,2)*1000,new_state(:,1)*1000)
% xlabel("downrange -- vbar [m]")
% ylabel("altitude -- rbar [m]")
% title("Hop 4 from -300 m vbar chaser location to hold 4 start location (-300 m --> -20 m)")

% find ending velocity of chaser
v_end_hop4_chaser_lvlh = new_state(end,4:6)'; % km/s
r_end_hop4_target_eci = new_state(7:9)';
v_end_hop4_target_eci = new_state(10:12)';

new_state_hop4 = new_state;

% ********************** PLOT IN ECI *****************************
% For all next ECI conversions
rho_eci = Q_eci_to_lvlh'*r_start_hop4_chaser_lvlh;
vrel_eci   = Q_eci_to_lvlh' * (v_start_hop4_chaser_lvlh + cross(omega_RSW, r_start_hop4_chaser_lvlh));% [km/s]
r_start_hop4_chaser_eci = r_initial_GEO + rho_eci;        % absolute chaser position in ECI [km]
v_start_hop4_chaser_eci = v_initial_GEO + vrel_eci;       % absolute chaser velocity in ECI [km/s]

timespan = [0 T_GEO]; %secs, multiplication to get it into seconds
state = [r_start_hop4_chaser_eci; v_start_hop4_chaser_eci];
[~,state_hop4_chaser_eci] = ode45(@twobodymotion,timespan,state,options,muearth);

% figure
% addEarth();
% hold on
% plot3(stateGEO_initial(:,1),stateGEO_initial(:,2),stateGEO_initial(:,3),"--","LineWidth",1.5);
% xlabel("x [km]")
% ylabel("y [km]")
% zlabel("z [km]")
% title("Hop 4 in ECI (-300 m --> -20 m)")
% plot3(state_hop4_chaser_eci(:,1),state_hop4_chaser_eci(:,2),state_hop4_chaser_eci(:,3),"LineWidth",1)
% plot3(r_initial_GEO(1,1), r_initial_GEO(2,1), r_initial_GEO(3,1), 'o', 'MarkerSize', 8, 'MarkerFaceColor', 'k')
% plot3(state_hop4_chaser_eci(1,1), state_hop4_chaser_eci(1,2), state_hop4_chaser_eci(1,3), '*', 'MarkerSize', 14, 'Color', 'm')
% plot3(state_hop4_chaser_eci(end,1), state_hop4_chaser_eci(end,2), state_hop4_chaser_eci(end,3), '*', 'MarkerSize', 14, 'Color', 'b')
% legend(" ","Target Trajectory","Chaser Trajectory","Starting Position of Target","Starting position of Chaser","Ending Position of Chaser")
% axis equal 
% grid on

%% Hold 4: Football

a = -0.02;
syms dvx
eq1 = a == 2*(dvx/n);
sol = solve(eq1,dvx);
dvx = double(sol);

r_start_hold4_chaser_lvlh = r_end_hop4_chaser_lvlh; %LVLH km
v_start_hold4_chaser_lvlh = [dvx;	0;	0]; % km/s

tspan = [0, T_GEO];
state = [r_initial_GEO_eci; v_initial_GEO_eci; r_start_hold4_chaser_lvlh; v_start_hold4_chaser_lvlh]; 
options = odeset('RelTol',1e-8,'AbsTol',1e-8);

[~,relativeOrbits] = ode45(@relativeMotion_Roshan,tspan,state,options,mu);

RC = [relativeOrbits(:,7),relativeOrbits(:,8),relativeOrbits(:,9)];
VC = [relativeOrbits(:,10),relativeOrbits(:,11),relativeOrbits(:,12)];

VCnorm = vecnorm(VC, 2, 2);

% figure('Name','Relative Distance')
% plot(RC(:,2)*1000,RC(:,1)*1000)
% hold on
% plot(0, 0, '*')
% grid on
% xlabel('Y Dist (m)')
% ylabel('X Dist (m)')

new_state_hold4 = new_state;

% ********************** PLOT IN ECI *****************************
% For all next ECI conversions
rho_eci = Q_eci_to_lvlh'*r_start_hold4_chaser_lvlh;
vrel_eci   = Q_eci_to_lvlh' * (v_start_hold4_chaser_lvlh + cross(omega_RSW, r_start_hold4_chaser_lvlh));% [km/s]
r_start_hold4_chaser_eci = r_initial_GEO + rho_eci;        % absolute chaser position in ECI [km]
v_start_hold4_chaser_eci = v_initial_GEO + vrel_eci;       % absolute chaser velocity in ECI [km/s]

timespan = [0 T_GEO]; %secs, multiplication to get it into seconds
state = [r_start_hold4_chaser_eci; v_start_hold4_chaser_eci];
[~,state_hold4_chaser_eci] = ode45(@twobodymotion,timespan,state,options,muearth);

% figure
% addEarth();
% hold on
% plot3(stateGEO_initial(:,1),stateGEO_initial(:,2),stateGEO_initial(:,3),"--","LineWidth",1.5);
% xlabel("x [km]")
% ylabel("y [km]")
% zlabel("z [km]")
% title("Hold 4 in ECI - Football")
% plot3(state_hold4_chaser_eci(:,1),state_hold4_chaser_eci(:,2),state_hold4_chaser_eci(:,3),"LineWidth",1)
% plot3(r_initial_GEO(1,1), r_initial_GEO(2,1), r_initial_GEO(3,1), 'o', 'MarkerSize', 8, 'MarkerFaceColor', 'k')
% plot3(state_hold4_chaser_eci(1,1), state_hold4_chaser_eci(1,2), state_hold4_chaser_eci(1,3), '*', 'MarkerSize', 14, 'Color', 'm')
% plot3(state_hold4_chaser_eci(end,1), state_hold4_chaser_eci(end,2), state_hold4_chaser_eci(end,3), '.', 'MarkerSize', 14, 'Color', 'g')
% legend(" ","Target Trajectory","Chaser Trajectory","Starting Position of Target","Starting position of Chaser","Ending Position of Chaser")
% axis equal 
% grid on

%% Maneuver 5: Final Approach (vbar)

% CHANGE T_GEO to whatever time we want to make exactly 10 days

% y0 = 20/1000; % km
vc = 20/1000/((10*24*60*60 - 10*T_GEO)); % km/s
time_final_approach = (10*24*60*60 - 10*T_GEO);
r_start_final_approach_chaser_lvlh = [0;-20/1000;0]; % km
v_start_final_approach_chaser_lvlh = [0;vc;0]; % km/s

% make sure we only have vbar rel velocity component

r_start_final_approach_target_eci = r_initial_GEO_eci;
v_start_final_approach_target_eci = v_initial_GEO_eci;

% time_final_approach = 1*24*60*60; % CHANGE LATER, whatever we need to get to exacly 10 days

% dydt = v_bar_approach(time, state, muearth, n)
% ODE initial state (12x1): [x y z xdot ydot zdot rtx rty rtz vtx vty vtz]
state =  [r_start_final_approach_chaser_lvlh;
          v_start_final_approach_chaser_lvlh;
          r_start_final_approach_target_eci;
          v_start_final_approach_target_eci];
[~, new_state] = ode45(@v_bar_approach,[0 time_final_approach],state,options,muearth,n);

% plot in LVLH (downrange vs altitude)
% figure
% plot(new_state(:,2)*10^3,new_state(:,1)*10^3)
% xlabel("downrange -- vbar [m]")
% ylabel("altitude -- rbar [m]")
% title("Final approach from chaser location to capture (-20 m --> 0 m)")

new_state_finalapproach = new_state;

% ********************** PLOT IN ECI *****************************
% For all next ECI conversions
rho_eci = Q_eci_to_lvlh'*r_start_final_approach_chaser_lvlh;
vrel_eci   = Q_eci_to_lvlh' * (v_start_final_approach_chaser_lvlh + cross(omega_RSW, r_start_final_approach_chaser_lvlh));% [km/s]
r_start_final_approach_chaser_eci = r_initial_GEO + rho_eci;        % absolute chaser position in ECI [km]
v_start_final_approach_chaser_eci = v_initial_GEO + vrel_eci;       % absolute chaser velocity in ECI [km/s]

timespan = [0 time_final_approach]; %secs, multiplication to get it into seconds
state = [r_start_final_approach_chaser_eci; v_start_final_approach_chaser_eci];
[~,state_fa_chaser_eci] = ode45(@twobodymotion,timespan,state,options,muearth);

stateGEO_fa = new_state(:,7:9);

figure
addEarth();
hold on
plot3(stateGEO_fa(:,1),stateGEO_fa(:,2),stateGEO_fa(:,3),"--","LineWidth",1.5);
xlabel("x [km]")
ylabel("y [km]")
zlabel("z [km]")
title("Final Approach in ECI (-20 m --> 0 m)")
plot3(state_fa_chaser_eci(:,1),state_fa_chaser_eci(:,2),state_fa_chaser_eci(:,3),"LineWidth",1)
plot3(stateGEO_fa(end,1), stateGEO_fa(end,2), stateGEO_fa(end,3), 'o', 'MarkerSize', 8, 'MarkerFaceColor', 'k')
plot3(state_fa_chaser_eci(1,1), state_fa_chaser_eci(1,2), state_fa_chaser_eci(1,3), '*', 'MarkerSize', 14, 'Color', 'm')
plot3(state_fa_chaser_eci(end,1), state_fa_chaser_eci(end,2), state_fa_chaser_eci(end,3), '*', 'MarkerSize', 14, 'Color', 'b')
legend(" ","Target Trajectory","Chaser Trajectory","Starting Position of Target","Starting position of Chaser","Ending Position of Chaser")
axis equal 
grid on
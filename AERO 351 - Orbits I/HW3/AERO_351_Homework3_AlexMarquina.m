
% AERO 351 HW 3
% Alex Marquina
% 10/25/24

close all;              % Clears all
clear;                  % Clears workspace
clc;                    % Clears command line

%% Constants
muEarth = 398600; % 
rEarth = 6378; % [km]
g = 9.807/1000; % [km/s^2]

%% Problem 4.15

% Givens
zp_415 = 300; % [km]
inc_415 = 35; % [degrees]
RAAN_415 = 130; % [degrees]
omega_415 = 115; % [degrees]
ecc_415 = 1.5; % [~]
theta_415 = 0; % [degrees] At perigee

% Calculate perigee radius
rp_415 = rEarth + zp_415; % [km]

% Call COEs to RV function
[R_GCF_415,V_GCF_415,R_Perifocal_415,V_Perifocal_415] = COEstoRV(muEarth,rp_415,ecc_415,inc_415,RAAN_415,omega_415,theta_415);

% Combine vectors
GCF_415 = [R_GCF_415 V_GCF_415]; % [km km/s]
Perifocal_415 = [R_Perifocal_415 V_Perifocal_415]; % [km km/s]

% Display results
disp('<<<<< Problem 4.15 >>>>>');
disp(' ');
disp('Perifocal radius vector [km]');
disp(R_Perifocal_415);
disp('Perifocal velocity vector [km/s]');
disp(V_Perifocal_415);
disp('GCF radius vector [km]');
disp(R_GCF_415);
disp('GCF velocity vector [km/s]');
disp(V_GCF_415);
disp(' ');

%% Problem 5.6

% Givens
r1_56 = [5644; 2830; 4170]; % [km]
r2_56 = [-2240; 7320; 4980]; % [km]
t_pass_56 = 20*60; % [s] Time between the r vectors

% Calculate inclination
perp_vec_56 = cross(r1_56,r2_56); % Perpendicular vector to the orbit plane
inc_56 = acosd(perp_vec_56(end)/norm(perp_vec_56)); % [degrees] Inclination


%% Problem 6.8

% Givens
z_1_68 = 300; % [km] Initial circular orbit
z_2_68 = 3000; % [km] Final circular orbit
r_1_68 = z_1_68 + rEarth; % [km]
r_2_68 = z_2_68 + rEarth; % [km]

% Calculate semimajor axis of the transfer orbit, energy of both  final orbits
epsilon_1_68 = -0.5*(muEarth/r_1_68); % [J/kg] Specific energy of the first orbit
v_1_68 = sqrt(muEarth/r_1_68); % [km/s] Velocity of the first orbit
epsilon_2_68 = -0.5*(muEarth/r_2_68); % [J/kg] Specific energy of the second orbit
v_2_68 = sqrt(muEarth/r_2_68); % [km/s] Velocity of the second orbit
a_68 = (r_1_68 + r_2_68)/2; % [km] Semimajor axis of the transfer orbit
ecc_68 = (r_2_68 - r_1_68)/(r_1_68 + r_2_68); % [~] Eccentricity of the transfer orbit

% Calculate angular momentum of the transfer orbit and velocities at perigee and apogee
h_68 = sqrt(muEarth*a_68*(1-ecc_68^2)); % [km^2/s^2]
vpt_68 = h_68/r_1_68; % [km/s] Speed at perigee of the transfer orbit
vat_68 = h_68/r_2_68; % [km/s] Speed at apogeee of the transfer orbit

% Calculate delta V
deltaV1 = vpt_68 - v_1_68; % [km/s]
deltaV2 = vat_68 - v_2_68; % [km/s]
deltaVTot = abs(deltaV2) + abs(deltaV1); % [km/s]

% Calculate time of the transfer
t_transfer_68 = pi*a_68^(3/2)/sqrt(muEarth); % [s]
t_transfer_68 = t_transfer_68/60; % [min]

% Display results
disp('<<<<< Problem 6.8 >>>>>');
disp(' ');
disp(['Total Delta V: ',num2str(deltaVTot),' km/s']);
disp(['Transfer time: ',num2str(t_transfer_68),' minutes']);
disp(' ');

%% Problem 6.23

% Givens
rp1_623 = 8100; % [km] Radius of perigee
ra1_623 = 18900; % [km] Radius of apogee
theta_intersect_623 = 45; % [degrees] Location of spacecraft B
theta1_623 = 150; % [degrees] Location of spacecraft C

% Calculate semimajor axis, eccentricity, angular momentum, and velocities at the point of intersection for the first orbit
a1_623 = (ra1_623 + rp1_623)/2; % [km] Semimajor axis
ecc1_623 = (ra1_623 - rp1_623)/(ra1_623 + rp1_623); % [~] Eccentricity
h1_623 = sqrt(muEarth*a1_623*(1-ecc1_623^2));
r1_623 = a1_623*(1-ecc1_623^2)/(1+ecc1_623*cosd(theta_intersect_623)); % [km] Radius of intersection
v1_azi_623 = muEarth*(1+ecc1_623*cosd(theta_intersect_623))/h1_623; % [km/s] Azimuthat velocity at B
v1_radial_623 = muEarth*ecc1_623*sind(theta_intersect_623)/h1_623; % [km/s] Radial velocity at B

% Calculate time since perigee for spacecraft B for the first orbit
[~,~,~,t1_623] = timeSincePerigeeEllipse(muEarth,theta_intersect_623,ecc1_623,a1_623);

% Calculate time since perigee for spacecraft C for the first orbit
[~,~,~,t2_623] = timeSincePerigeeEllipse(muEarth,theta1_623,ecc1_623,a1_623);

% Calculate orbit period and time since perigee of the first orbit at S/C C
T1_623 = 2*pi*a1_623^(3/2)/sqrt(muEarth); % [s] Period of the first orbit

% Calculate the required period of orbit 2, semimajor axis, energy, and velocity at S/C B
T2_623 = T1_623 + t1_623 - t2_623; % [s] Period requirement of the second orbit
a2_263 = (T2_623*sqrt(muEarth)/(2*pi))^(2/3); % [km] Semimajor axis of the second orbit
epsilon2_623 = -0.5*muEarth/a2_263; % [J/kg] Energy of the second orbit
v_2_623 = sqrt((epsilon2_623 + muEarth/r1_623)*2); % [km/s] Velocity of the first orbit at intersection

% Calculate delta v requirement for the phasing maneuver
deltaVTransfer_623 = 2*sqrt((v1_azi_623 - v_2_623)^2 + v1_radial_623^2); % [km/s]

% Display results
disp('<<<<< Problem 6.23 >>>>>');
disp(' ');
disp(['Delta V required: ',num2str(deltaVTransfer_623),' km/s']);
disp(' ');

%% Problem 6.25

% Givens
zp_625 = 1270; % [km] Perigee altitude
rp_625 = zp_625 + rEarth; % [km] Radius of perigee
vp_625 = 9; % [km/s] Perigee speed
ecc2_625 = 0.4; % [~] Target eccentricity
theta_625 = 100; % [degrees] True anomaly location for maneuver
thetaEnd_625 = 180; % [degrees] Theta B (NH with common apse)

% Calculate angular momentum, ecc, and radius at 100 degrees true anomaly for the first orbit
h1_625 = vp_625*rp_625; % [km^2/s^2]
ecc1_625 = (h1_625^2)/(muEarth*rp_625) - 1; % [~] Eccentricity of the first orbit
r1_625 = (h1_625^2)/(muEarth*(1 + ecc1_625*cosd(theta_625))); % [km]
v1_azi_625 = (muEarth/h1_625)*(1+ecc1_625*cosd(theta_625)); % [km/s] Azimuthal velocity of the first orbit at 100 degrees
v1_radial_625 = (muEarth/h1_625)*ecc1_625*sind(theta_625); % [km/s]
v1_625 = sqrt(v1_radial_625^2 + v1_azi_625^2); % [km/s]
gamma1_625 = atand(v1_radial_625/v1_azi_625); % [degrees] Flight path angle at 100 degrees

% Calculate radius of apogee and semimajor axis for the final orbit
a2_625 = rp_625/(1-ecc2_625); % [km] Semimajor axis of the final orbit
ra2_625 = a2_625*(1 + ecc2_625); % [km] Radius of apogee for the final orbit

% Determine requirements for the transfer orbit
h3_625 = sqrt(muEarth*r1_625*(1 + ecc2_625*cosd(theta_625)));
v2_radial_625 = (muEarth/h3_625)*(ecc2_625)*sind(theta_625);
v2_azi_625 = (muEarth/h3_625)*(1 + (ecc2_625*cosd(theta_625)));
gamma2_625 = atand(v2_radial_625/v2_azi_625);
v2_625 = sqrt(v2_azi_625^2 + v2_radial_625^2);

% Calculate delta V and gamma
deltaGamma = gamma2_625 - gamma1_625; % [degrees]
deltaV_625 = sqrt(v1_625^2 + v2_625^2 - 2*v1_625*v2_625*cosd(deltaGamma));

% Display results
disp('<<<<< Problem 6.25 >>>>>');
disp(' ');
disp(['Delta gamma: ',num2str(deltaGamma),' degrees']);
disp(['Delta V: ',num2str(deltaV_625),' km/s']);
disp(' ');


%% Problem 6.44

% Givens
z1_644 = 300; % [km] Starting altitude
z2_644 = 600; % [km] Target altitude
inc2_644 = 20; % [degrees] Target inclination
r1_644 = z1_644 + rEarth; % [km] Initial radius
r2_644 = z2_644 + rEarth; % [km] Target radius

% Find velocity of the first and second orbit
v1_644 = sqrt(muEarth/r1_644); % [km/s] Velocity of the first orbit
v2_644 = sqrt(muEarth/r2_644); % [km/s] Velocity of the second orbit

% Find transfer Hohmann orbit semimajor axis, eccentricity, angular momentum, and velocities at apogee and perigee 
a3_644 = (r1_644 + r2_644)/2; % [km] Semimajor axis of the transfer orbit
ecc3_644 = (r2_644 - r1_644)/(r1_644 + r2_644); % [~] Eccentricity of the transfer orbit
h3_644 = sqrt(muEarth*a3_644*(1-ecc3_644^2)); % [km^2/s^2] Angular momentum of the tramsfer orbit
vp3_644 = h3_644/r1_644; % [km/s] Perigee velocity of the transfer orbit
va3_644 = h3_644/r2_644; % [km/s] Apogee velocity of the transfer orbit

%%%%%% Part a %%%%%

deltaV_a_HT = (vp3_644 - v1_644) + (v2_644 - va3_644); % [km/s] Hohmann delta V requirement
deltaV_a_inc = 2*v2_644*sind(inc2_644/2); % [km/s] Inclination delta V requirement
deltaV_a_tot = deltaV_a_inc + deltaV_a_HT; % [km/s] Total delta V requirement

%%%%% Part b %%%%%

deltaV_b_HT = (vp3_644 - v1_644); % [km/s]
deltaV_b_inc = sqrt(va3_644^2 + v2_644^2 - 2*va3_644*v2_644*cosd(inc2_644)); % [km/s]
deltaV_b_tot = deltaV_b_inc + deltaV_b_HT; % [km/s]

%%%%% Part c %%%%%

deltaV_c_inc = sqrt(vp3_644^2 + v1_644^2 - 2*vp3_644*v1_644*cosd(inc2_644)); % [km/s]
deltaV_c_HT = v2_644 - va3_644; % [km/s]
deltaV_c_tot = deltaV_c_HT + deltaV_c_inc; % [km/s]

% Display results
disp('<<<<< Problem 6.44 >>>>>');
disp(' ');
disp(['Part a delta V: ',num2str(deltaV_a_tot),' km/s']);
disp(['Part b delta V: ',num2str(deltaV_b_tot),' km/s']);
disp(['Part c delta V: ',num2str(deltaV_c_tot),' km/s']);
disp(' ');

%% Problem 6.47

% Givens
r1_647 = [436; 6083; 2529]; % [km]
v1_647 = [-7.34; -0.5125; 2.497]; % [km/s]
coastTime = [0 89*60]; % [s] Coast time for the spacecraft
Isp = 300; % [s] Thruster Isp
m = 1000; % [kg] Satellite mass
F = 10; % [kN] Thruster force
burnTime = [0, 120]; % [s] Engine burn time

% Propagate the coasting orbit
state1 = [r1_647; v1_647];
options = odeset('RelTol',1e-8,'AbsTol',1e-8);
[timenew1,statenew1] = ode45(@twobodymotion,coastTime,state1,options,muEarth);

% Propagate the continuous thrust function
state2 = [statenew1(end,1); statenew1(end,2); statenew1(end,3); statenew1(end,4); statenew1(end,5); statenew1(end,6); m];
[timenew2,statenew2] = ode45(@continuousThrust,burnTime,state2,options,muEarth,F,Isp,g);

R_647 = [statenew2(end,1); statenew2(end,2); statenew2(end,3)];
V_647 = [statenew2(end,4); statenew2(end,5); statenew2(end,6)];

[a_647, ecc_647, TrueAnom_647, inc_647, RAAN_647, AoP_647, h_647, epsilon_647, ra_647, rp_647, T_647, E, F, Me, Mh, t] = RVtoCOEs(R_647, V_647, muEarth);

% Calculate the time to hit max altitude and the max altitude
zmax_647 = ra_647 - rEarth; % [km] Maximum altitude
ra_time = T_647/2;
total_time = (ra_time - t) + 120 + 89*60;

% Display results
disp('<<<<< Problem 6.47 >>>>>');
disp(' ');
disp(['Max altitude: ',num2str(zmax_647),' km']);
disp(['Total time: ',num2str(total_time/60/60),' hours']);
disp(' ');


%% Functions

% COEs to RV function
function [R,V,rPerifocal,vPerifocal,h] = COEstoRV(mu,rp,ecc,inc,RAAN,AoP,theta)
    h = sqrt(mu*(1+ecc)*rp);
    rPerifocal = (h^2/mu)/(1+ecc*cosd(theta)) *[cosd(theta);sind(theta);0];
    vPerifocal = (mu/h)*[-sind(theta);ecc+cosd(theta);0];

    Cz_RAAN = [cosd(RAAN)  sind(RAAN) 0
               -sind(RAAN) cosd(RAAN) 0
                    0           0     1];

    Cx_Inc = [1     0         0
              0 cosd(inc) sind(inc)
              0 -sind(inc) cosd(inc)];

    Cz_AoP = [cosd(AoP)   sind(AoP) 0
               -sind(AoP) cosd(AoP) 0
                   0          0     1];

    DCM = Cz_AoP*Cx_Inc*Cz_RAAN;
    
    % GCF frame
    R = DCM'*rPerifocal;
    V = DCM'*vPerifocal;

end

% Time since periapse passage
function [E,Me,n,t] = timeSincePerigeeEllipse(mu,theta,ecc,a)
    E = 2*atan(sqrt((1-ecc)/(1+ecc))*tand(theta/2)); % [rad] Eccentric anomaly
    Me = E - ecc*sin(E); % [rad] Mean anomaly
    n = sqrt(mu/(a^3)); % [rad/s] Mean motion
    t = Me/n; % [s] Time since perigee
end

% ODE 45 function for coast
function dstate = twobodymotion(~,state,muearth)

    % Define variables
    x = state(1);
    y = state(2);
    z = state(3);
    dx = state(4);
    dy = state(5);
    dz = state(6);

    % Magnitude of the position vector
    r = norm([x y z]);
    
    ddx = -muearth*x/r^3;
    ddy = -muearth*y/r^3;
    ddz = -muearth*z/r^3;

    dstate = [dx; dy; dz; ddx; ddy; ddz];

end

% ODE 45 function for continuous thrust
function dstate = continuousThrust(~,state,mu,T,Isp,g)
    
    % Define variables
    x = state(1);
    y = state(2);
    z = state(3);
    dx = state(4);
    dy = state(5);
    dz = state(6);
    m = state(7);

    % Equations of motion and mass
    r = norm([x y z]);
    v = norm([dx dy dz]);
    ddx = (-mu*x)/r^3 + (T*dx/(m*v));
    ddy = (-mu*y)/r^3 + (T*dy/(m*v));
    ddz = (-mu*z)/r^3 + (T*dz/(m*v));
    dm = -T/(Isp*g);
    
    dstate = [dx; dy; dz; ddx; ddy; ddz; dm];
end



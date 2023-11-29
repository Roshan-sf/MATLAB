%Roshan Jaiswal-Ferri
%Aero 215 HW4 Orbit Transfer: 11/28/23

%%

close all;      %Clears all
clear all;      %Clears Workspace
clc;            %Clears Command Window

%% Defining Variables

R = [6161.56, 454.07, -2537.72]; %km
V = [0.376, 7.391, 2.224]; %km/s
mu = 398600; %in km^3/S^2

%% Part 1: Calculating COEs for satelite

[a,e,nu,i,RAAN,w,p] = hw3_orbitalCOEs_Jaiswal_ferriRoshan(R,V,mu);

% Converting Rad to Deg
nu = rad2deg(nu);
RAAN = rad2deg(RAAN);
i = rad2deg(i);
w = rad2deg(w);
p = p/3600; %seconds to hours

% Displaying Results
disp('Part 1 - Results for Initial Orbit:')
disp(['     Eccentricity: ', num2str(e), ' unitless'])
disp(['     Semi-Major Axis: ', num2str(a), ' km'])
disp(['     True Anomaly: ', num2str(nu), ' deg'])
disp(['     Inclination: ', num2str(i), ' deg'])
disp(['     RAAN: ', num2str(RAAN), ' deg'])
disp(['     Argument of Periapsis: ', num2str(w), ' deg'])
disp(['     Period: ', num2str(p), ' hours'])
disp(' ')
disp(' ')

%% Part 2: Geostationary Orbit

%New Vectors:
R_geo = [42157, 0, 0]; %km
V_geo = [0, 3.07, 0]; %km/s

[a1,e1,nu1,i1,RAAN1,w1,p1] = hw3_orbitalCOEs_Jaiswal_ferriRoshan(R_geo,V_geo,mu);

% Converting Rad to Deg
nu1 = rad2deg(nu1);
RAAN1 = rad2deg(RAAN1);
i1 = rad2deg(i1);
w1 = rad2deg(w1);
p1 = p1/3600; %seconds to hours

% Displaying Results
disp('Part 2 - Results for Final Geostationary Orbit:')
disp(['     Eccentricity: ', num2str(e1), ' unitless'])
disp(['     Semi-Major Axis: ', num2str(a1), ' km'])
disp(['     True Anomaly: ', num2str(nu1), ' deg'])
disp(['     Inclination: ', num2str(i1), ' deg'])
disp(['     RAAN: ', num2str(RAAN1), ' deg'])
disp(['     Argument of Periapsis: ', num2str(w1), ' deg'])
disp(['     Period: ', num2str(p1), ' hours'])
disp(' ')
disp(' ')

%% Part 3: Four Burn Transfer
%sqrt(mu/norm(R)) = velocity to be at for a circular orbit R circular =
%a*(1+e);

epsilon = ((norm(V)^2))/(2) - ((mu)/(norm(R))); %specific mech energy

%Burn 1: Circularizing orbit (e of 0)

Rcircular = a*(1+e); %r circular which points to apoapsis
Vgiven = sqrt(2*(mu/norm(Rcircular)+epsilon)); 
Vcircular = sqrt(mu/norm(Rcircular)); %Velocity final circular orbit
deltaVCirc = abs(Vgiven - Vcircular); %delta V circular

%Burn 2: Changing inclination to 0

deltaTheta = i - 0; %in degrees!!
deltaVinc = 2*Vcircular*sind(deltaTheta/2); %delta v inclination change sind cuz degrees!!

%Burn 3: Hohman burn 1

ah1 = ((norm(Rcircular)+norm(R_geo))/2); %Calculating semi major axis for first transfer
eph1 = -((mu)/(2*ah1));
Vh1 = sqrt(2*((mu/norm(Rcircular)+eph1)));
deltaVh1 = abs(Vh1 - Vcircular);

%Burn 4: Hohman burn 2

eph2 = -((mu)/(2*ah1));
Vh2 = sqrt(2*((mu/norm(R_geo)+eph2)));
deltaVh2 = abs(Vh2 - norm(V_geo));

%Displaying Delta V results

disp('Part 3 - Delta V for Transfer Burns:')
disp(['     Burn 1 (Circularization): ', num2str(deltaVCirc), ' km/s']); %km/s
disp(['     Burn 2 (Inclination): ', num2str(deltaVinc), ' km/s']); %km/s
disp(['     Burn 3 (Hohmann 1): ', num2str(deltaVh1), ' km/s']); %km/s
disp(['     Burn 4 (Hohmann 2): ', num2str(deltaVh2), ' km/s']); %km/s
disp(['     Total Delta V: ', num2str(deltaVCirc + deltaVinc + deltaVh1 + deltaVh2), ' km/s'])
disp(' ')
disp(' ')

%% Part 4: Three Burn Transfer

%Burn 1: Circularizing orbit (e of 0)

Rcircular = a*(1+e); %r circular which points to apoapsis
Vgiven = sqrt(2*(mu/norm(Rcircular)+epsilon)); 
Vcircular = sqrt(mu/norm(Rcircular)); %Velocity final circular orbit
deltaVCirc = abs(Vgiven - Vcircular); %delta V circular

%Burn 2: Hohmann Burn 1

ah1 = ((norm(Rcircular)+norm(R_geo))/2); %Calculating semi major axis for first transfer
eph1 = -((mu)/(2*ah1));
Vh1 = sqrt(2*((mu/norm(Rcircular)+eph1)));
deltaVh1 = abs(Vh1 - Vcircular);

%Burn 3: Combined Plane Change (Hohmann Burn 2 Inclination)

Vcpc = sqrt((Vh2^2) + (norm(V_geo)^2) - 2*(Vh2)*(norm(V_geo))*cosd(deltaTheta));

%Display Results

disp('Part 4 - Delta V for Transfer Burns with CPC:')
disp(['     Burn 1 (Circulation): ', num2str(deltaVCirc), ' km/s']); %km/s
disp(['     Burn 2 (Hohmann 1): ', num2str(deltaVh1), ' km/s']); %km/s
disp(['     Burn 3 (Combined Plane Change): ', num2str(Vcpc), ' km/s']); %km/s
disp(['     Total Delta V: ', num2str(deltaVCirc + deltaVh1 + Vcpc), ' km/s'])




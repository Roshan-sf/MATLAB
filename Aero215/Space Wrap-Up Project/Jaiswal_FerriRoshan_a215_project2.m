%Roshan Jaiswal-Ferri
%Aero 215 Space Wrap-Up Project: 12/07/23

%%

close all;      %Clears all
clear all;      %Clears Workspace
clc;            %Clears Command Window

%% Defining Variables

CD.R = [-406.663, -4186.877, -5059.146]; %km
CD.V = [7.386, -2.178, 1.1889]; %km/s

ISS.R = [5648.682, -2337.321, 2943.766]; %km
ISS.V = [-0.208, 5.799, 5.008]; %km/s

mu = 398600; %in km^3/S^2

%% Part 1: Calculating COEs & Comparing Apogee, Perigee, and Inclination

% Calculating COEs

%For Cargo Dragon
[CD.a,CD.e,CD.nu,CD.i,CD.RAAN,CD.w,CD.p] = hw3_orbitalCOEs_Jaiswal_ferriRoshan(CD.R,CD.V,mu);
%For ISS
[ISS.a,ISS.e,ISS.nu,ISS.i,ISS.RAAN,ISS.w,ISS.p] = hw3_orbitalCOEs_Jaiswal_ferriRoshan(ISS.R,ISS.V,mu);

% Calculating Apoapsis and Periapsis
% Apogee = a*(1+e), Perigee =  a*(1-e) 

CD.ap = (CD.a)*(1+CD.e); %km
CD.pe = (CD.a)*(1-CD.e); %km

ISS.ap = (ISS.a)*(1+ISS.e); %km
ISS.pe = (ISS.a)*(1-ISS.e); %km

% Displaying Results
disp('Part 1 - Apogee, Perigee, & Inclination:')
disp(' ')
disp(' Results for ISS: ')
disp(['     Apogee: ', num2str(ISS.ap), ' km'])
disp(['     Perigee: ', num2str(ISS.pe), ' km'])
disp(['     Inclination: ', num2str(rad2deg(ISS.i)), ' deg'])
disp(' Results for CD: ')
disp(['     Apogee: ', num2str(CD.ap), ' km'])
disp(['     Perigee: ', num2str(CD.pe), ' km'])
disp(['     Inclination: ', num2str(rad2deg(CD.i)), ' deg'])
disp(' ')
disp(' ')


%The Perigee of the ISS and Cargo Dragon are very similar to eachother as
%well as the Apogee. They are almost identical, this is
%probably because they are supposed to rendezvous together. The orbit of
%the Cargo Dragon was specifically altered to match the ISS to meet.

%It will probably rendezvous closer perigee point as it is closer between
%each orbit than the apogee.

%The inclination is almost exactly the same, also probably because they are
%supposed to dock to eachother. The Cargo Dragons inclination is off
%because its you cannot have a perfect launch.

%% Part 2: Circularizing Cargo Dragons orbit

CD.se = ((norm(CD.V)^2))/(2) - ((mu)/(norm(CD.R))); %specific mech energy

CD.Vap = sqrt(2*(mu/norm(CD.ap)+CD.se)); %Velocity at apogee
CD.Vc = sqrt(mu/norm(CD.ap)); %Velocity after circular change
CD.dVc = abs(CD.Vap - CD.Vc); %delta V circular km/s

disp('Part 2 - Delta V for Circulization:')
disp([' Delta V For Circulization: ', num2str(CD.dVc), ' km/s'])
disp(' ')
disp(' ')

%% Part 3: 

g0 = 9.80665; %m/s^2
Isp = 316; %seconds
Mf = 12568; %Mass final in kg
CD.dVc = CD.dVc*1000; %convert from km/s to m/s

Mi = Mf*exp((CD.dVc)/(Isp*g0));
Mp = abs(Mf - Mi); %kg

disp('Part 3 - Mass of CD and Propellant:')
disp([' Mass Initial: ', num2str(Mi), ' kg'])
disp([' Mass Final: ', num2str(Mf), ' kg'])
disp([' Mass Propellant: ', num2str(Mp), ' kg'])

MpV = Mp/1000;

disp([' Propellant Volume: ', num2str(MpV), ' m^3'])









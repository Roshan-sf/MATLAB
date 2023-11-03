function [a,e,nu,i,RAAN,w] = hw3_orbitalCOEs_Jaiswal_ferriRoshan(R,V,mu)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

%% Part 1: Defining Variables

RM = norm(R); %Magnitude of R
VM = norm(V); %Magnitude of V

ui = [1,0,0];
uj = [0,1,0];
uk = [0,0,1];
h = cross(R,V);

uiM = norm(ui); %the magnitudes of the values above
ujM = norm(uj);
ukM = norm(uk);
hM = norm(h);


%% Part 2: Initial Calculations for lalter

ep = ((VM^2)/2)-((mu)/RM); %Calculating Epsilon (specific mechanical energy) in J/kg


%% Part 3: Calculating semi-major axis

a = -((mu)/(2*ep)); %in km

%% Genreal equation calculation for period

p = (2*pi)*sqrt((a^3)/(mu)); %period of orbit in seconds

%% Part 4: Calculating eccentricity
eV = (1/mu)*((((VM^2)-((mu)/(RM)))*R)-(dot(R,V)*V));

e = norm(eV);

%% Part 5:

i = acos((dot(uk,h))/((hM)*(ukM))); %in rad not deg


%% Part 6: RAAN

n = cross(uk,h); %projection of momentum vector in orbital plane
nM = norm(n);

RAAN = acos((dot(ui,n))/((uiM)*(nM)));

%% Part 7: Argument of Periapsis

w = acos(dot(n,eV)/(nM*e));

%% Part 8: nu (or theta)

nu = acos(dot(eV,R)/(e*RM));


































end
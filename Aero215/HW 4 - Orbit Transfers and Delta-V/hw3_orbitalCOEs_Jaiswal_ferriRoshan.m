function [a,e,nu,i,RAAN,w,p] = hw3_orbitalCOEs_Jaiswal_ferriRoshan(R,V,mu)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

%% Part 1: Defining Variables

RM = norm(R); %Magnitude of R
VM = norm(V); %Magnitude of V

ui = [1,0,0];
uj = [0,1,0];
uk = [0,0,1];
h = cross(R,V);
h2 = dot(R,V);

uiM = norm(ui); %the magnitudes of the values above
ujM = norm(uj);
ukM = norm(uk);
hM = norm(h);


%% Part 2: Initial Calculations for later

ep = ((VM^2)/2)-((mu)/RM); %Calculating Epsilon (specific mechanical energy) in J/kg


%% Part 3: Calculating semi-major axis

a = -((mu)/(2*ep)); %in km

%% Genreal equation calculation for period

p = (2*pi)*sqrt((a^3)/(mu)); %period of orbit in seconds

%% Part 4: Calculating eccentricity
eV = (1/mu)*((((VM^2)-((mu)/(RM)))*R)-(dot(R,V)*V)); %eccentricity vector is from origin to point of periapsis 

e = norm(eV);

%% Part 5: inclination in rad

i = acos((dot(uk,h))/((hM)*(ukM))); %in rad not deg


%% Part 6: RAAN in rad

n = cross(uk,h); %projection of momentum vector in orbital plane and node line?
nM = norm(n);

if n(2) >= 0
    RAAN = acos((dot(ui,n))/((uiM)*(nM))); %original equation
else
    RAAN = (2*pi)-(acos((dot(ui,n))/((uiM)*(nM))));
end

%% Part 7: Argument of Periapsis in rad

if eV(3) >= 0 %k component of eccentricity vector (height)
    w = acos(dot(n,eV)/(nM*e));
else
    w = (2*pi)-(acos(dot(n,eV)/(nM*e)));
end

%% Part 8: nu (or theta) true anomaly in rad

if h2 >= 0 %dot product of R and V idk what it represents
    nu = acos(dot(eV,R)/(e*RM));
else
    nu = (2*pi)-(acos(dot(eV,R)/(e*RM)));
end




end
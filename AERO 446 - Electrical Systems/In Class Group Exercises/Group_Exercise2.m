%% Roshan Jaiswal-Ferri, Stefan Rosu, Jack Schafer, Jordan Powell
%Section - 01
%Aero 446 GE2: 4/22/25

%% Workspace Prep

%warning off
format long     %Allows for more accurate decimals
close all;      %Clears all
clear all;      %Clears Workspace
clc;            %Clears Command Window

%% S/C & Env Data

mu = 398600;
Re = 6378;
altVec = 300:20:2000;
Se = 1367;
lty = 10; %lifetime of 10 yrs

Rsc = Re + altVec;
beta = 0;
P = 2*pi*sqrt(Rsc.^3/mu);

tPayload = 2*60;
tComms = 5*60;

power.Payload = 138;
power.Structure = 30;
power.Thermal = 80;
power.Power = 25;
power.Comms = 40;
power.Computer = 25;
power.const = power.Structure + power.Thermal + power.Power + power.Computer;

%% Question 1

fe = (1/pi) * asin( sqrt((Re./Rsc).^2 - sin(beta)^2) / cos(beta) );

tEclipse = P.*fe;

figure('Name','Eclipse Time vs Orbital Altitude')
plot(altVec, tEclipse)
ylabel('Eclipse Time (s)')
xlabel('Altitude (km)')
title('Eclipse Time vs Orbital Altitude')

%% Question 2

tDay = P - tEclipse;

figure('Name','Day Time vs Orbital Altitude')
plot(altVec, tDay)
ylabel('Day Time (s)')
xlabel('Altitude (km)')
title('Day Time vs Orbital Altitude')

%% Question 3

powGenTime = P./2;

figure('Name','Power Gen Time vs Orbital Altitude')
plot(altVec, powGenTime)
ylabel('Power Gen Time (s)')
xlabel('Altitude (km)')
title('Power Gen Time vs Orbital Altitude')

%% Question 4

E_const_night = power.const .* tEclipse;

E_comms = power.Comms * tComms;

% Total energy at night (J)
E_night = E_const_night + E_comms;

% Average discharge power over eclipse
power.night = E_night ./ tEclipse;

E_const_day = power.const .* tDay;

E_payload = power.Payload * tPayload;

E_charge = (power.const.*((P./2)-tEclipse)) + E_night;

% Total energy needed during the day (J)
E_day = E_const_day + E_payload + E_charge;

% Average power needed over the daylight duration
power.day = E_day ./ tDay;

power.charge = (E_charge)/(P./2);


figure('Name','Power Day')
plot(altVec, power.day)
ylabel('Power Req Day (s)')
xlabel('Altitude (km)')
title('Power Req Day vs Orbital Altitude')

%% Question 5

figure('Name','Power Req Night vs Orbital Altitude')
plot(altVec, power.night)
ylabel('Power Req Night (s)')
xlabel('Altitude (km)')
title('Power Req Night vs Orbital Altitude')

%% Question 6

n = 1; %Assume 100% eff
Fp = 1;
theta = 45; %degrees

arrayA = power.day ./ (Se * n * Fp * cosd(theta));

%% Question 7

Dyr = 0.002; %degredation per year in percent
lifedegredation = (1-Dyr)^lty;
peol = power.day.*lifedegredation; %in W
Af = power.day./peol;
Adiff = Af.*arrayA;

figure('Name','Solar Array Area vs Altitude')
plot(altVec, arrayA)
hold on
grid on
plot(altVec,Adiff)
xlabel('Altitude (km)')
ylabel('Required Solar Array Area (m^2)')
legend('BOL','EOL',Location='best')
title('Solar Array Area vs Altitude')

%% Question 8

n_bat = 0.8;
DOD = 0.5;

battCap = (power.night.*(tEclipse))/(DOD*n_bat); %Watt seconds
battCap = battCap./3600; %watt hours

figure('Name','Battery Capacity vs Altitude')
plot(altVec,battCap)
ylabel('Battery Capacity (Wh)')
xlabel('Altitude (km)')
title('Battery Capacity vs Altitude')

%% Question 9

cost = missionCost(altVec, arrayA);

figure('Name','Cost vs Altitude')
plot(altVec,cost)
ylabel('Cost ($)')
xlabel('Altitude (km)')
title('Cost vs Altitude')

%% Question 10

[~, pos] = min(cost);
disp(['Lowest mission cost of [$', num2str(cost(pos)), '], at altitude of: ', num2str(altVec(pos)), 'km'])

%% Question 11

sunAngle = linspace(pi/2,(3*pi)/2, 86); %rads
%power.gen = power.day ./ (Se * n * Fp .* -cos(sunAngle));

zero = zeros(1,86);
space1 = linspace(0,pi/2, 86); 
Psp2 = -arrayA(1)*(Se*n*Fp*cos(sunAngle));
space2 = linspace((3*pi)/2, 2*pi, 86);

Psp = [zero,Psp2,zero];
location = [space1,sunAngle,space2];

figure('Name','Power Gen on Orbit')
plot(location,Psp)
xticks(0:pi/2:2*pi)  % Set tick locations
xticklabels({'0','\pi/2','\pi','3\pi/2','2\pi'})  % Set tick labels
xlabel('Angle (rad)')
ylabel('Power Generation (Watt)')
xlabel('Radial Position (radians)')
title('Power Gen on Orbit')

%% Functions

function Cost = missionCost(altVec, arrayA)
    a = 1.2;
    b = 1.5;
    k_altitude = 50000;
    k_solar = 5000;
    C_baseline = 3e6;

    Cost = C_baseline + k_altitude * (altVec .^ a) + k_solar * (arrayA .^ b);
end
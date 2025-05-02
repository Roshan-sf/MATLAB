%% Roshan Jaiswal-Ferri
%Section - 01
%Aero 356 Midterm 1: 4/23/25

%% Workspace Prep

%warning off
format long     %Allows for more accurate decimals
close all;      %Clears all
clear all;      %Clears Workspace
clc;            %Clears Command 

%% Problem 1

%Orbital Parameters
Re = 6378.18; %km
mu = 398600;
alt = 500; %km
inc = 23.5; %deg
rp = alt + Re;
ra = rp;
Se = 1367; %w/m^2
sb = 5.67*10^-8;

%mission parameters
lty = 7;
lt = 7*365.25; %7 years
Fp = 0.85;
sf = 0.87; %average sunlight fraction (percent time in sun each orbit)
a = (rp + ra)/2;
p = (2*pi)*sqrt((a^3)/(mu)); %period in seconds

% Question 1
opd = (24*3600)/p; %orbits per day
pm = p/60; %period in minutes %it is 94.6 minutes, good for LEO

avg_st = pm*sf; %average total time spent within sunlight (sun time)
et = pm - avg_st; %average time in eclipse per orbit
etd = et*opd; %average time in eclipse spent per day in minutes

disp('Average time spent in eclipse per day (minutes): ')
disp(num2str(etd))
disp(' ')

% Question 2
no = 2.660; %nominal operating voltage for Voc
BOL = -5.9*10^-3; %milli-voltage drop per degree C
EOL = -6.5*10^-3;

syms x x2
eq1 = BOL/no == x/100;
eq2 = EOL/no == x2/100;

soln = solve([eq1,eq2],x,x2);

degB = double(soln.x);
degE = double(soln.x2);

disp('Percent Voc Drop per deg C at BOL:')
disp(num2str(degB))
disp('Percent Voc Drop per deg C at EOL:')
disp(num2str(degE))
disp('Percent degredation difference between BOL-EOL')
disp(num2str(abs(degE-degB)))
disp(' ')

% Question 3

%solar array properties
cella = 26.6; %area per cell in cm^2
array = 60*3; %cells per array
arrayAcm = (array*cella);
arrayA = arrayAcm/10000; %solar array area m^2 !!front or back only!!
abs = 0.92;
E = 0.85;

Pabs = abs*Se*arrayA;

syms T
eq3 = Pabs == sb*E*(arrayA)*T^4;
soln2 = solve(eq3,T);
Temp = double(soln2(2));
TempC = Temp-273.15;

disp('Equilibrium Temp in C:')
disp(num2str(TempC))
disp(' ')

% Question 4

std_power = 135.3; %std power generated per area (mW/cm^2) @ 28C
Tempdiff = TempC - 28;
degbol = (Tempdiff*-degB)/100; %drop in performance at current temp
%pbol = arrayAcm*std_power*degbol*Fp; 

n0 = 28.1/100;
n = n0*(1+(degB/100)*(Tempdiff)); %the bug in the code making too much power is here :/
pbol = Se*n*Fp*arrayA; %power generated in mW at BOL per 1 panel

disp('Power generated in W at BOL:')
disp(num2str(2*(pbol)))
disp(' ')

% Question 5

Dyr = 0.005; %degredation per year in percent
lifedegredation = (1-Dyr)^lty;
peol = pbol*lifedegredation; %in mW

disp('Power generated in W at EOL:')
disp(num2str(2*(peol)))
disp(' ')

syms F %size factor
eq4 = 215 == (peol*2)*F;
soln3 = solve(eq4,F);
sizef = double(soln3);

%find new area
newArraycm = arrayAcm*sizef;
cells = (newArraycm/cella)-(2*60*3);

disp('Minimum whole new cells needed to reach 215W at EOL & equi temp:')
disp(num2str(ceil(cells)))
disp(' ')

%% Problem 2

% Question 1

syms r
circ = 30200 == 2*pi*r;
soln4 = solve(circ,r);
newR = double(soln4);

disp('New Earth Radius (km):')
disp(num2str(newR))
disp(' ')

syms a2
daysec = 86400;
period = daysec == (2*pi)*sqrt((a2^3)/(mu)); %period in seconds
a2 = solve(period,a2);
r = double(a2(1))/2;

disp('CEO altitude (km):')
disp(num2str(r))
disp(' ')

% Question 2

FOR = 2*asind(newR/r);

disp('Field of Regard for s/c in CEO orbit (deg):')
disp(num2str(FOR))







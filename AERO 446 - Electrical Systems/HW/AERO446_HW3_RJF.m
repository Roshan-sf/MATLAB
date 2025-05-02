%% Roshan Jaiswal-Ferri
%Section - 01
%Aero 446 HW3: 4/25/25

%% Workspace Prep

%warning off
format long     %Allows for more accurate decimals
close all;      %Clears all
clear all;      %Clears Workspace
clc;            %Clears Command Window

%% Vars

%Subsystem 	Power (W)
Payload = 138;
Structure = 3;
Thermal = 10;
Power = 5;
TTC = 36; %used for up and downlink
OBC = 5; %On Board Computer
ADCS = 30;

contPower = Structure + Thermal + Power + OBC;

%orbit parameters
p = 90*60; %period in seconds
dT = 67; %day temp in C
eT = -65; %Eclipse temp in C
dl = 1/3; %amt of time spent in eclipse
dli = 2/3;
Xd = 0.8;
Xe = 0.9;
Se = 1366; %w/m^2

%solar cell
cellArea = 0.06*0.05; %m
powDens = 135.3; %mw/cm^2
n0 = 0.302; %std efficiency at 28 C
Voc = 0.99; %volts
Tcoeff = -6.7*10^-3; % V/deg C
Fp = 1; %assume packing fraction of 1 since unlisted

%% Problem 1

%% Question 1.1

Tc = (Tcoeff/Voc);

disp(['Solar Cell Temp coefficient: ', num2str(Tc*100), ' %/C'])
disp(' ')

%% Question 1.2

Td = p*dli;
Te = p*dl;

%account for eclipse
Je1 = contPower*Te; %energy in joules
Je2 = TTC*3*60; %energy in joules for 3 min downlink
Je3 = ADCS*30; %joules for 30sec ADCS activation night
Je = Je1 + Je2 + Je3;

%account for daylight
Jd1 = contPower*Td;
Jd2 = Payload*2*60; %joules for payload for 2 min per orbit
Jd3 = ADCS*30; %joules for 30sec ADCS activation day
Jd = Jd1 + Jd2 + Jd3;

Po = orbitPower(Je,Jd,Te,Td,Xe,Xd); %power for each orbit

disp(['Avg s/c power required per orbit: ', num2str(Po), ' W'])
disp(' ')
%% Question 1.3

n = n0*(1+Tc*(dT-28)); %using day temps
Jd = Se*n*Fp*cellArea*Td; %power in joules
Pe = Jd/p;

disp(['Avg power generated per orbit per cell: ', num2str(Pe), ' W'])
disp(' ')

%% Question 1.4

cellCount = ceil(Po/Pe);
disp(['Cells needed to supply power: ', num2str(cellCount)])
disp(' ')

%% Question 1.5

arrayA = cellCount*cellArea/2;
disp(['Area required per panel: ', num2str(arrayA), ' m^2'])
disp(' ')

%% Question 1.6

% These solar cells are adequate for a box-panel type s/c because the
% panels are relatively small while generating enough power to stay
% opertional. They also function properly within the required temperature
% range.

%% Question 1.7

% Potential improvements to the power system of this spacecraft could
% include powering the ADCS for less time or using one that uses less power
% over a similar activation time. Other suggestions would be to increase 
% data storage on the space craft to allow for less frequent downlinks.

%% Problem 2

%% Qustion 2.1

% The recent and sad death of the opportunity mars rover was because of a
% dust storm covering the solar panels in dust and lowering power input so
% much that it ran out of battery. By switching to a RTG NASA eliminated a
% variable that was out of their control and directly addressed the 
% previous point of failure. Now they can control their power input and do
% not have to rely on the sun which could be covered. RTGs will also 
% increase the lifespan of the rover, which has been shown to be pretty 
% long. 

%% Question 2.2

% You can use capacitors or batteries to handle the peak loads of a system
% while making sure to size the solar array for appropriate power
% generation.

%% Question 2.3

% Peak to peak power means energy is stored during low demand and matched
% at high demand. This is useful because it allows you to store extra power
% and can improve the efficiency of the system.

%% Question 2.4

% Solar panels will have better performance coming out of an eclipse
% because their operating temperature will genreally be cooler and closer
% to the optimal performance temp. The temperature differential slowly
% rises over the day time of the orbit and decreases efficiency.

%% Question 2.5

% Solar panel degredation is definitely worse in LEO because of earth's
% atmosphere-components like a lot of AO, outgassing/particulates in the 
% higher atmosphere. They will also go through many more thermal cycles in
% LEO which can cause fatigue.

%% Problem 3

% 1) False
% 2) False
% 3) True
% 4) False
% 5) True

%% Functions

function [power] = orbitPower(Je,Jd,Te,Td,Xe,Xd)
    c1 = Je/Xe;
    c2 = Jd/Xd;
    power = (c1 + c2) / Td;
end


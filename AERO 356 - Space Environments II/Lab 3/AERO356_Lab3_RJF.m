%% Roshan Jaiswal-Ferri
%Section - 01
%Aero 356 Lab 3: 5/22/25

%% Workspace Prep

warning off
format long     %Allows for more accurate decimals
close all;      %Clears all
clear all;      %Clears Workspace
clc;            %Clears Command 

%% Reading Data

POLYData = readtable("Data\LangmuirData.csv");
Vp = POLYData.V;
Ip = -((POLYData.V_I_MV).*(1e-3))./220;
min2 = min(Ip);
Ipfix = Ip - min2;

UCLAData = readtable("Data\LangmuirProfessional.xlsx");
V = UCLAData.voltage;
I = UCLAData.current;
uplift = abs(min(I)*1.00001);
Ifix = I + uplift;

%% Finding Floating Potential

% Sort the absolute values and get their indices
[~, idx] = sort(abs(I));

% Take the indices of the 10 values closest to zero
numClosest = min(8, length(idx));  % Ensure we don't exceed vector length
closestIdx = idx(1:numClosest);

% Select those values from V
Vfv = V(closestIdx);

Vf = abs(mean(Vfv));
disp(['Floating Potential Voltage: ', num2str(Vf)])
disp(' ')

%% Ion Saturation Current

[~, idx] = min(abs(V+5)); %it was chosen that the saturation zone ended ~5V

IionSat = I(1:idx);

IonSat = abs(mean(IionSat));
disp(['Ion Saturation Voltage: ', num2str(IonSat)])
disp(' ')

%% Finding Plasma Temp

%11V to 28V for slope
[~, idx_startSlope] = min(abs(V-15));
[~, idx_endSlope] = min(abs(V-29));
% logData = log10(Ifix);

Islope = I(idx_startSlope:idx_endSlope);
Vslope = V(idx_startSlope:idx_endSlope);
p = polyfit(Vslope,log(Islope),1);

Te = p(1)^-1; %eV

%% Find Voltage Potential

[~, idx_startElectronSat] = min(abs(V-40));

IElectronSat = I(idx_startElectronSat:end);
VElectronSat = V(idx_startElectronSat:end);
p2 = polyfit(VElectronSat,log(IElectronSat),1);

m = p(1);
b = p(2);

m2 = p2(1);
b2 = p2(2);

syms x
eq = (m*x) + b == (m2*x) + b2;
soln = solve(eq,x);
VpUCLA = double(soln); %V

%% Plasma Density

D = 0.0508/100;
r = 0.5*D;
L = 0.1905/100;
Ap = 2*pi*r*L; %area m^2
q = 1.60217663e-19; %coulumbs
k = 1.380649e-23;
m = 1.6726e-27; %kg of protons in xenon
mp = m*(131.293); %amount of protons and nuetrons

syms n
eq2 = IonSat == 0.6*q*n*sqrt((k*Te*11606)/mp)*Ap;
soln2 = solve(eq2,n);
np = double(soln2);

%% Collected Data

% Sort the absolute values and get their indices
[~, idx2] = sort(abs(Ip));

% Take the indices of the 10 values closest to zero
numClosest2 = min(10, length(idx2)); %Talk about 20 value, how many min points, half the data (only 60 pts)
closestIdx2 = idx2(1:numClosest2); 

% Select those values from V
Vfv2 = Vp(closestIdx2);

Vf2 = abs(mean(Vfv2));
disp(['Floating Potential Voltage: ', num2str(Vf2)])
disp(' ')

%% Ion Saturation Current

% this is saying -5 volts, its really neg but it has absolute

[~, idx2] = min(abs(Vp+5)); %it was chosen that the saturation zone ended ~5V

IionSat2 = Ip(1:idx2);

IonSat2 = abs(mean(IionSat2));
disp(['Ion Saturation Voltage: ', num2str(IonSat2)])
disp(' ')

%% Finding Plasma Temp

%40V to 51V for slope
[~, idx_startSlope2] = min(abs(Vp-40));
[~, idx_endSlope2] = min(abs(Vp-51));

Islope2 = Ip(idx_startSlope2:idx_endSlope2);
Vslope2 = Vp(idx_startSlope2:idx_endSlope2);
p3 = polyfit(Vslope2,log(Islope2),1);

Te = p3(1)^-1 %eV

%% Find Voltage Potential

[~, idx_startElectronSat2] = min(abs(Vp-70));
IElectronSat2 = Ip(idx_startElectronSat2:end);
VElectronSat2 = Vp(idx_startElectronSat2:end);
p2 = polyfit(VElectronSat2,log(IElectronSat2),1);
m = p(1);
b = p(2);
m2 = p2(1);
b2 = p2(2);
syms x
eq = (m*x) + b == (m2*x) + b2;
soln = solve(eq,x);
Vp = double(soln);


%%

% figure()
% plot(V,I)
% grid on
% 
% figure()
% plot(Vp,Ip)
% grid on
% 
% figure()
% semilogy(Vp,Ipfix)
% grid on





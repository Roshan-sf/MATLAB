%% Roshan Jaiswal-Ferri
%Section - 01
%Aero 356 Lab 3: 5/22/25

%% Workspace Prep

%warning off
format long     %Allows for more accurate decimals
close all;      %Clears all
clear all;      %Clears Workspace
clc;            %Clears Command 

%% Reading Data

data = readtable("Data\LangmuirProfessional.xlsx");
V = data.voltage;
I = data.current;
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
Vp = double(soln); %V

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

%%

% figure()
% plot(V,I)
% grid on
% 
% figure()
% plot(V,Ifix)
% grid on
% 
% figure()
% semilogy(V,Ifix)
% grid on




%% Roshan Jaiswal-Ferri
%Section - 01 
%Aero 303 Thermocouple Plotting: 3/10/25

%% Workspace Prep

format long     %Allows for more accurate decimals
close all;      %Clears all
clear all;      %Clears Workspace
clc;            %Clears Command Window

%% 

Al = readtable("Lab3_1A/aluminum_temperature.dat");
Cu = readtable("Lab3_1A/copper_temperature.dat");
St = readtable("Lab3_1A/steel_temperature.dat");
HotC = readtable("Lab3_1A/hot_calibration.dat");
ColdC = readtable("Lab3_1A/ice_calibration.dat");

%%

figure("Name",'Cold Calibration')
plot(ColdC.Var2,ColdC.Var4)
xlabel('Time (s)')
ylabel('Temperature (F)')
title('Cold Calibration (Ice)')

figure("Name",'Hot Calibration')
plot(HotC.Var2,HotC.Var4)
xlabel('Time (s)')
ylabel('Temperature (F)')
title('Hot Calibration (Boiling Water)')


figure("Name",'Overlayed Calibration Data')
plot(HotC.Var2,HotC.Var4)
hold on
grid on
plot(ColdC.Var2,ColdC.Var4)
xlabel('Time (s)')
ylabel('Temperature (F)')
title('Overlayed Calibration Data')
legend('Hot Calibration','Cold Calibration')


%% Error

avgEC = abs(mean(ColdC.Var4) - 32);
avgEH = abs(mean(HotC.Var4) - 212);
avgE = mean([avgEH,avgEC]); %temp avg

avgEC = abs(mean(ColdC.Var4) - 32)/32;
avgEH = abs(mean(HotC.Var4) - 212)/212;
avgEP = mean([avgEH,avgEC]); %percent avg

disp(['Avg Percent Error: ', num2str(avgEP)]);
disp(['Avg Degree Error (F): ', num2str(avgE)]);

%% Heat Plot in Rod

figure("Name",'Rod Temperature vs Time')
plot(Al.Var2,Al.Var4)
hold on
grid on
plot(Cu.Var2,Cu.Var4)
plot(St.Var2,St.Var4)
xlabel('Time (s)')
ylabel('Temperature (F)')
title('Rod Temperature vs Time')
legend('Aluminum Rod','Copper Rod','Steel Rod',Location='best')

%% 1D Heat transf - Convection—ROD to AIR
%Cu:
Ts = 322.594; %hottest temp of copper rod in K (121 F)
Tinf = 327.594; %Temp of steam in K (130 F)
Tf = (Ts+Tinf)/2;
B = 1/Tf;
v = 1.5049e-5; %interpolated and found from table 28
g = 9.8;
Lc = 0.1778; %meters
Pr = 0.71;
Gr = (g*B*(Tinf-Ts)*Lc^3)/(v^2);
Nu = 0.59*(Gr*Pr)^(0.25);
he = Lc;
r = 0.003175;
A = 2*pi*r*he;
kair = 0.024594444444444443;
h = (Nu*kair)/Lc;
qCu = h*A*(Tinf-Ts);

%Al
Ts = 308.20555555555552019; %hottest temp of al rod in K (F)95 
Tinf = 327.594; %Temp of steam in K (130 F)
Tf = (Ts+Tinf)/2;
B = 1/Tf;
v = 1.5049e-5; %interpolated and found from table 28
g = 9.8;
Lc = 0.1778; %meters
Pr = 0.71;
Gr = (g*B*(Tinf-Ts)*Lc^3)/(v^2);
Nu = 0.59*(Gr*Pr)^(0.25);
he = Lc;
r = 0.003175;
A = 2*pi*r*he;
kair = 0.024594444444444443;
h = (Nu*kair)/Lc;
qAl = h*A*(Tinf-Ts);

%St
Ts = 318.03889; %hottest temp of al rod in K (112.8 F)
Tinf = 327.594; %Temp of steam in K (130 F)
Tf = (Ts+Tinf)/2;
B = 1/Tf;
v = 1.5049e-5; %interpolated and found from table 28
g = 9.8;
Lc = 0.1778; %meters
Pr = 0.71;
Gr = (g*B*(Tinf-Ts)*Lc^3)/(v^2);
Nu = 0.59*(Gr*Pr)^(0.25);
he = Lc;
r = 0.003175;
A = 2*pi*r*he;
kair = 0.024594444444444443;
h = (Nu*kair)/Lc;
qSt = h*A*(Tinf-Ts);

%% 1D Heat transf - Conduction—WATER to ROD
%Cu


%% Workspace Prep

format long     %Allows for more accurate decimals
close all;      %Clears all
clear all;      %Clears Workspace
clc;            %Clears Command Window

Pr = 0.7496;
T0 = 1916.667;
rho = 0.268; %kg/m^3
gamma = 1.25;
R = 287;
T = T0/(1+((gamma-1)/2));
a = sqrt(gamma*R*T);
v = a;
D = 0.003175; %1/8 in in meters
L = D;
mu = 48.445/10^6;
Re = (rho*v*L)/mu;

Nu = 0.023*(Re^0.8)*(Pr^0.4);


x = 1;

















% k_al = 237;  %Thermal conductivity of Aluminum (W/mK)
% k_cu = 401;  %Thermal conductivity of Copper (W/mK)
% k_cs = 50.2; %Thermal conductivity of Carbon Steel (W/mK)
% 
% L = 12.25/39.37;  %Rod length (5 inches = 0.127 m)
% A = pi * (0.00635)^2;  % Cross-sectional area (for 1/4 inch rod)
% T_base = 212;  % Base temperature in Fahrenheit
% T_air = 55;  % Ambient air temperature in Fahrenheit
% 
% % Convert temperatures to Celsius
% T_base_C = (T_base - 32) * 5/9;
% T_air_C = (T_air - 32) * 5/9;
% 
% % Heat Transfer Rate Calculation using Fourier's Law
% Q_al = (k_al * A / L) * (T_base_C - T_air_C);
% Q_cu = (k_cu * A / L) * (T_base_C - T_air_C);
% Q_cs = (k_cs * A / L) * (T_base_C - T_air_C);
% 
% % Display results
% disp(['Heat Transfer Rate for Aluminum Rod: ', num2str(Q_al), ' W']);
% disp(['Heat Transfer Rate for Copper Rod: ', num2str(Q_cu), ' W']);
% disp(['Heat Transfer Rate for Carbon Steel Rod: ', num2str(Q_cs), ' W']);

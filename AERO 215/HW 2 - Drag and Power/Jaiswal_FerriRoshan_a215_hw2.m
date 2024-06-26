%Roshan Jaiswal-Ferri
%Aero 215 HW 2: 10/19/23

%%

close all;      %Clears all
clear all;      %Clears Workspace
clc;            %Clears Command Window

%% Part 1 Defining Variables

%Variables of aircraft parameters
aircraft.WS = 8.27; %in lbf/ft^2
aircraft.AR = 5.63; %unitless
aircraft.CD0 = .0275; %unitless
aircraft.e = .78; %unitless
aircraft.CLM = 1.4; %unitless
aircraft.S = 127; %ft^2
aircraft.rho =  .0023769; %rhoSealevel in slugs/ft^3


aircraft_7500 = aircraft;
aircraft_7500.rho = 0.00189722581; %rho7500ft in slugs/ft^3

aircraft_landing = aircraft;
aircraft_landing.CD0 = .0975; %landing config Drag coef
aircraft_landing.CLM = 1.7; %landing config max lift

aircraft_heavy = aircraft;
aircraft_heavy.WS = 10.39; %heavy wing loading in lbf/ft^2

%% Part 2 Confirming Hand Calculations
disp('%% Part 2:')
disp(' ')
v_max = 70; %in knots 

%Set 1 of calculations at sea level

[V, D, P] = hw2_dragPower_Jaiswal_FerriRoshan(v_max, aircraft); %Baseline sea level

V1 = V(end);
D1 = D(end);
P1 = P(end);
disp('  %% Part 2, Set 1:')
disp(['     Velocity (kts): ',num2str(V1),', Drag Force (lbf): ',num2str(D1),','])
disp(['     Drag Power (hp): ',num2str(P1)])
disp(' ')

%Set 2 of calculations at 7,500ft

[V_7500, D_7500, P_7500] = hw2_dragPower_Jaiswal_FerriRoshan(v_max, aircraft_7500); %Baseline 7500ft

V1 = V_7500(end);
D1 = D_7500(end);
P1 = P_7500(end);
disp('  %% Part 2, Set 2:')
disp(['     Velocity (kts): ',num2str(V1),', Drag Force (lbf): ',num2str(D1),','])
disp(['     Drag Power (hp): ',num2str(P1)])
disp(' ')

%Set 3 of calculations at landing config

[V_landing, D_landing, P_landing] = hw2_dragPower_Jaiswal_FerriRoshan(v_max, aircraft_landing); %landing config sea level

V1 = V_landing(end);
D1 = D_landing(end);
P1 = P_landing(end);
disp('  %% Part 2, Set 3:')
disp(['     Velocity (kts): ',num2str(V1),', Drag Force (lbf): ',num2str(D1),','])
disp(['     Drag Power (hp): ',num2str(P1)])
disp(' ')

%Set 4 of calculations at heavy wing loading config

[V_heavy, D_heavy, P_heavy] = hw2_dragPower_Jaiswal_FerriRoshan(v_max, aircraft_heavy); %heacy loading sea level

V1 = V_heavy(end);
D1 = D_heavy(end);
P1 = P_heavy(end);
disp('  %% Part 2, Set 4:')
disp(['     Velocity (kts): ',num2str(V1),', Drag Force (lbf): ',num2str(D1),','])
disp(['     Drag Power (hp): ',num2str(P1)])
disp('     Yay they match!')
disp(' ')




%% Part 3 Plotting

v_max = 140; %in knots

[V, D, P] = hw2_dragPower_Jaiswal_FerriRoshan(v_max, aircraft); %Baseline sea level
[V_7500, D_7500, P_7500] = hw2_dragPower_Jaiswal_FerriRoshan(v_max, aircraft_7500); %Baseline 7500ft
[V_landing, D_landing, P_landing] = hw2_dragPower_Jaiswal_FerriRoshan(v_max, aircraft_landing); %landing config sea level
[V_heavy, D_heavy, P_heavy] = hw2_dragPower_Jaiswal_FerriRoshan(v_max, aircraft_heavy); %heacy loading sea level

figure; %creating a figure with three graphs of Temp, Pres, and Density vs altitude
subplot(1,2,1)
plot(V, D,'red') 
hold on;
plot(V_7500, D_7500,'c')
plot(V_landing, D_landing,'y')
plot(V_heavy, D_heavy,'g')
xlabel('Velocity (kts)');
ylabel('Drag Force (lbf)')
title('Drag Force vs Velocity')
grid on;
legend('Baseline', '7500ft', 'Landing','Heavy', location = 'best')

subplot(1,2,2)
plot(V, P,'red') 
hold on;
plot(V_7500, P_7500,'c')
plot(V_landing, P_landing,'y')
plot(V_heavy, P_heavy,'g')
xlabel('Velocity (kts)');
ylabel('Drag Power (hp)')
title('Drag Power vs Velocity')
grid on;
legend('Baseline', '7500ft', 'Landing','Heavy', location = 'best')

%% Part 4 Energy Height

%Variables:
q = 0.3048; %Conversion ratio for f to m
e = 3.28084; %Conversion ratio for m/s to ft/s
ma = .7; %mach speed
g = 32.1741; %g in ft/s
R = 287.058; %gas constant 
gma = 1.4; %Gamma isentropic gas parameter
Eh = 0; %Energy height in ft
h0 = 0; %altitude at sea level (in either ft or m)
hft = 10000; %Height in ft
h = hft*q; %converting ft to m for stadm function


[T, ~, ~] = stdatm_Jaiswsal_FerriRoshan(h); %calling the stdatm function for T
sps = sqrt(gma*T*R); %speed of sound in m/s
spsft = sps*e; %speed of sound in ft/s

[T, ~, ~] = stdatm_Jaiswsal_FerriRoshan(h0); %calling the stdatm function for T
spss = sqrt(gma*T*R); %speed of sound at sea level in m/s
spssft = spss*e; %speed of sound in ft/s


v = spsft*ma; % calculating speed in ft/s
Ze = ((.5)*((v^2)/(g)))+hft; %in ft
vf = sqrt((2*Ze*g)-(2*h0*g)); %calculating v at sea level in ft/s
maf = vf/spssft;
disp('%% Part 4:')
disp(['Speed of Sound at 10,000 ft in ft/s: ',num2str(spsft)])
disp('Mach number of aircraft at 10,000 ft: 0.7,');
disp(['Speed of Sound at sea level in ft/s: ',num2str(spssft)])
disp(['Mach number of aircraft at sea level: ',num2str(maf)])
disp(['Calculated energy height in ft: ',num2str(Ze)])










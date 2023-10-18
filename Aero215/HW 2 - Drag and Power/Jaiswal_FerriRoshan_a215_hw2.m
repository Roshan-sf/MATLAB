%Roshan Jaiswal-Ferri
%Aero 215 HW 2: 10/19/23

%%

close all;      %Clears all
clear all;      %Clears Workspace
clc;            %Clears Command Window

%% Part 1 Defining

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
aircraft_heavy.WS = 12.39; %heavy wing loading in lbf/ft^2

%% Part 2 Confirming Hand Calculations

v_max = 70; %in knots
[V, D, P] = hw2_dragPower_Jaiswal_FerriRoshan(v_max, aircraft); %Baseline sea level

V1 = V(end);
D1 = D(end);
P1 = P(end);

disp(['Velocity (kts): ',num2str(V1),', Drag Force (lbf): ',num2str(D1),', Drag Power (hp): ',num2str(P1)])
disp('Yay they match!')

%% Part 3 Plotting

v_max = 140; %in knots

[V, D, P] = hw2_dragPower_Jaiswal_FerriRoshan(v_max, aircraft); %Baseline sea level
[V_7500, D_7500, P_7500] = hw2_dragPower_Jaiswal_FerriRoshan(v_max, aircraft_7500); %Baseline sea level
[V_landing, D_landing, P_landing] = hw2_dragPower_Jaiswal_FerriRoshan(v_max, aircraft_landing); %Baseline sea level
[V_heavy, D_heavy, P_heavy] = hw2_dragPower_Jaiswal_FerriRoshan(v_max, aircraft_heavy); %Baseline sea level

figure; %creating a figure with three graphs of Temp, Pres, and Density vs altitude
subplot(1,2,1)
plot(V, D,'red') 
hold on;
plot(V_7500, D_7500,'c')
plot(V_landing, D_landing,'y')
plot(V_heavy, D_heavy,'g')
xlabel('Velocity (kts)');
ylabel('Drag Force (lbf)')
title('Velocity vs Drag Force')
grid on;
legend('Baseline', '7500ft', 'Landing','Heavy')

subplot(1,2,2)
plot(V, P,'red') 
hold on;
plot(V_7500, P_7500,'c')
plot(V_landing, P_landing,'y')
plot(V_heavy, P_heavy,'g')
xlabel('Velocity (kts)');
ylabel('Drag Power (hp)')
title('Velocty vs Drag Power')
grid on;
legend('Baseline', '7500ft', 'Landing','Heavy')

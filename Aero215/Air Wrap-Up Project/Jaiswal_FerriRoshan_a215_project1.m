%Roshan Jaiswal-Ferri
%Aero 215 Air Wrap-Up Project: 10/26/23

%%

close all;      %Clears all
clear all;      %Clears Workspace
clc;            %Clears Command Window

%% Part 1 Defining Variables

h = 0; %altitude in m

%Variables of aircraft parameters
RV10.EW = 1605; %aircraft empty weight in lbf
RV10.PS = 4; %max passengers
RV10.MF = 360; %Max fuel load in lbf
RV10.CA = 8000; %Cruise altitude in ft
RV10.CD0 = .0265; %unitless
RV10.AR = 6.81; %unitless
RV10.e = .79; %unitless
RV10.CLM = 1.15; %unitless
RV10.S = 148; %ft^2
RV10.VM = 175; %v max in kts

RV12.EW = 740; %aircraft empty weight in lbf
RV12.PS = 2; %max passengers
RV12.MF = 119; %Max fuel load in lbf
RV12.CA = 7500; %Cruise altitude in ft
RV12.CD0 = .0275; %unitless
RV12.AR = 5.63; %unitless
RV12.e = .78; %unitless
RV12.CLM = 1.4; %unitless
RV12.S = 127; %ft^2
RV12.VM = 125; %v max in kts

%CP is common parameters

CP.PSFC = .07; %Power Specific Fuel Consumption in gallons per hour per horsepower
CP.G = 5.04; %avg cost per gallon of fuel in USD
CP.GD = 6.01; %avg gas density in lbf per gallon
CP.PW = 170; %avg weight per passenger in lbf

%% Part 1: Calculating minimum drag lbf and its speed for each fully loaded aircraft at take off

[~, ~, rho] = stdatm_Jaiswsal_FerriRoshan(h); 

%Max Weight Results

[V, D, ~] = hw2_dragPower_Jaiswal_FerriRoshan(CP, RV10, rho);

[V2, D2, ~] = hw2_dragPower_Jaiswal_FerriRoshan(CP, RV12, rho);


%Min Weight Results

RV10.MF = 0; %fuel load in lbf
RV12.MF = 0; %fuel load in lbf

[V3, D3, ~] = hw2_dragPower_Jaiswal_FerriRoshan(CP, RV10, rho);

[V4, D4, ~] = hw2_dragPower_Jaiswal_FerriRoshan(CP, RV12, rho);

%Plotting Results

figure;
subplot(1,2,1)
plot(V, D,'red') 
hold on;
plot(V2, D2,'y')
xlabel('Velocity (kts)');
ylabel('Drag Force (lbf)')
title('Drag Force vs Velocity @Max Weight')
grid on;
legend('RV10', 'RV12', location = 'best')

subplot(1,2,2)
plot(V3, D3,'red') 
hold on;
plot(V4, D4,'y')
xlabel('Velocity (kts)');
ylabel('Drag Force (lbf)')
title('Drag Force vs Velocity @Min Weight')
grid on;
legend('RV10', 'RV12', location = 'best')

%Find Min Drag for each plane in each scenario 

%This is setting variable m1 to the minumum value in each drag vector
m1 = min(D); %Minimum drag at max weight  rv10
m2 = min(D2); %Minimum drag at max weight  rv12
m3 = min(D3); %Minimum drag at min weight  rv10
m4 = min(D4); %Minimum drag at min weight  rv12
%This finds the position value in each drag vector of the its minumum drag value
p = find(D==m1); 
p2 = find(D2==m2);
p3 = find(D3==m3);
p4 = find(D4==m4);
%This sets the velocity to the same position of the minimum drag in the velocity vector
v = V(p); %speed at minimum drag for max weight rv10
v2 = V2(p2); %speed at minimum drag for max weight rv12
v3 = V3(p3); %speed at minimum drag for min weight rv10
v4 = V4(p4); %speed at minimum drag for min weight rv12
%Display them
disp('Minimum Weight:')
disp(['  Minimum Drag For RV10: ', num2str(m3), ' lbf @', num2str(v3),' kts'])
disp(['  Stall Speed For RV10: ', num2str(V3(1)), ' kts'])
disp(['  Minimum Drag For RV12: ', num2str(m4), ' lbf @', num2str(v4),' kts'])
disp(['  Stall Speed For RV12: ', num2str(V4(1)), ' kts'])
disp(' ')
disp('Maximum Weight:')
disp(['  Minimum Drag For RV10: ', num2str(m1), ' lbf @', num2str(v),' kts'])
disp(['  Stall Speed For RV10: ', num2str(V(1)), ' kts'])
disp(['  Minimum Drag For RV12: ', num2str(m2), ' lbf @', num2str(v2),' kts'])
disp(['  Stall Speed For RV12: ', num2str(V2(1)), ' kts'])
disp(' ')
disp(' ')

%Making a table (this sucks)

%T = table('Size',[4 3],'VariableTypes',{'double','double','double'}, 'VariableNames',{'Optimal V (kts)', 'Min Drag (lbf)', 'Stall Speed (kts)'},'RowNames', {'RV10 Min Weight', 'RV12 Min Weight', 'RV10 Max Weight', 'RV20 Max Weight'});

Aircraft = {'RV10 Min Weight';'RV12 Min Weight';'RV10 Max Weight';'RV12 Max Weight'};
Optimal_Velocity_kts = [v3;v4;v;v2];
Drag_Force_lbf = [m3;m4;m1;m2];
Stall_Speed_kts = [V3(1);V4(1);V(1);V2(1)];

T = table(Aircraft,Optimal_Velocity_kts,Drag_Force_lbf,Stall_Speed_kts);

disp(T)


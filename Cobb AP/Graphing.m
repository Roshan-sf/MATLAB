%% Roshan Jaiswal-Ferri
% 
%2006 WRX STI Cobb AP Datalogs: 2/15/25

%% Workspace Prep

format long     %Allows for more accurate decimals
close all;      %Clears all
clear all;      %Clears Workspace
clc;            %Clears Command Window

%% Reading Data

city = readtable("datalog1.csv");
sport = readtable("datalog2.csv");

%% Assigning Vars

timeC = city.Time_sec_;
boostC = city.Boost_psi_;
clutchC = city.ClutchSw_on_off_;
coolantC = city.CoolantTemp_F_;
DAMC = city.DynAdvMult_DAM_;
feedbackknC = city.FeedbackKnock___;
EconomyC = city.FuelEconomy_mpg_;
GearC = city.GearPosition_gear_;
intakeTC = city.IntakeTemp_F_;
MAFC = city.MAF_g_s_;
manifoldPresC = city.ManAbsPress_psi_;
RPMC = city.RPM_RPM_;
throttlePosC = city.ThrottlePos___;
speedC = city.VehicleSpeed_mph_;

timeS = sport.Time_sec_;
boostS = sport.Boost_psi_;
clutchS = sport.ClutchSw_on_off_;
coolantS = sport.CoolantTemp_F_;
DAMS = sport.DynAdvMult_DAM_;
feedbackknS = sport.FeedbackKnock___;
EconomyS = sport.FuelEconomy_mpg_;
GearS = sport.GearPosition_gear_;
intakeTS = sport.IntakeTemp_F_;
MAFS = sport.MAF_g_s_;
manifoldPresS = sport.ManAbsPress_psi_;
RPMS = sport.RPM_RPM_;
throttlePosS = sport.ThrottlePos___;
speedS = sport.VehicleSpeed_mph_;

%% Plotting For City

% figure('Name','City Driving Data 1')
% plot(timeC,boostC,'b')
% hold on
% grid on
% plot(timeC,RPMC,'r')
% legend('Boost', 'RPM', Location='best')

figure('Name','City Driving Data 1')
yyaxis left  % Left y-axis for Boost
plot(timeC, boostC, 'b')
ylabel('Boost (psi)')
yyaxis right % Right y-axis for RPM
plot(timeC, RPMC, 'r')
ylabel('RPM')
grid on
legend('Boost', 'RPM', Location='best')
xlabel('Time (sec)')
title('City Driving Data 1')


figure('Name','City Driving Data 2')
plot(timeC,throttlePosC,'y')
hold on
grid on
plot(timeC,EconomyC,'c')
legend('throttle pos', 'fuel econ', Location='best')

% figure('Name','City Driving Data 3')
% plot(timeC,DAMC,'g')

%% Plotting For Sport (beginning w/ launch control)

figure('Name','Sport Driving Data 1')
yyaxis left  % Left y-axis for Boost
plot(timeS, boostS, 'b')
ylabel('Boost (psi)')
yyaxis right % Right y-axis for RPM
plot(timeS, RPMS, 'r')
ylabel('RPM')
grid on
legend('Boost', 'RPM', Location='best')
xlabel('Time (sec)')
title('Sport Driving Data 1')


figure('Name','Sport Driving Data 2')
plot(timeS,throttlePosS,'y')
hold on
grid on
plot(timeS,EconomyS,'c')
legend('throttle pos', 'fuel econ', Location='best')

figure('Name','City Driving Data 3')
plot(timeS,DAMS,'g')
legend('DAM', Location='best')

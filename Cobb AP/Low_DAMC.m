%% Roshan Jaiswal-Ferri
% 
%2006 WRX STI Cobb AP Datalogs: 2/15/25

%% Workspace Prep

format long     %Allows for more accurate decimals
close all;      %Clears all
clear all;      %Clears Workspace
clc;            %Clears Command Window
warning('off');

%% Reading Data

control1    = readtable("Datalogs\Low DAM\datalog.csv");
control2    = readtable("Datalogs\Low DAM\datalog2.csv");
reflash1    = readtable("Datalogs\Low DAM\datalog3.csv");
reflash1pt2 = readtable("Datalogs\Low DAM\datalog4.csv");
reflash2    = readtable("Datalogs\Low DAM\datalog5.csv");

%% Assigning Vars

% Control 1
timeC         = control1.Time_sec_;
boostC        = control1.Boost_psi_;
clutchC       = control1.ClutchSw_on_off_;
coolantC      = control1.CoolantTemp_F_;
DAMC          = control1.DynAdvMult_DAM_;
feedbackknC   = -1 .* control1.FeedbackKnock___;
EconomyC      = control1.FuelEconomy_mpg_;
GearC         = control1.GearPosition_gear_;
intakeTC      = control1.IntakeTemp_F_;
MAFC          = control1.MAF_g_s_;
manifoldPresC = control1.ManAbsPress_psi_;
RPMC          = control1.RPM_RPM_;
throttlePosC  = control1.ThrottlePos___;
speedC        = control1.VehicleSpeed_mph_;
wasteC        = control1.WastegateDuty___;

% Control 2
timeC2         = control2.Time_sec_;
boostC2        = control2.Boost_psi_;
clutchC2       = control2.ClutchSw_on_off_;
coolantC2      = control2.CoolantTemp_F_;
DAMC2          = control2.DynAdvMult_DAM_;
feedbackknC2   = -1 .* control2.FeedbackKnock___;
EconomyC2      = control2.FuelEconomy_mpg_;
GearC2         = control2.GearPosition_gear_;
intakeTC2      = control2.IntakeTemp_F_;
MAFC2          = control2.MAF_g_s_;
manifoldPresC2 = control2.ManAbsPress_psi_;
RPMC2          = control2.RPM_RPM_;
throttlePosC2  = control2.ThrottlePos___;
speedC2        = control2.VehicleSpeed_mph_;
wasteC2        = control2.WastegateDuty___;

% Reflash 1
timeR         = reflash1.Time_sec_;
boostR        = reflash1.Boost_psi_;
clutchR       = reflash1.ClutchSw_on_off_;
coolantR      = reflash1.CoolantTemp_F_;
DAMR          = reflash1.DynAdvMult_DAM_;
feedbackknR   = -1 .* reflash1.FeedbackKnock___;
EconomyR      = reflash1.FuelEconomy_mpg_;
GearR         = reflash1.GearPosition_gear_;
intakeTR      = reflash1.IntakeTemp_F_;
MAFR          = reflash1.MAF_g_s_;
manifoldPresR = reflash1.ManAbsPress_psi_;
RPMR          = reflash1.RPM_RPM_;
throttlePosR  = reflash1.ThrottlePos___;
speedR        = reflash1.VehicleSpeed_mph_;
wasteR        = reflash1.WastegateDuty___;

% Reflash 1 pt 2
timeR1pt2         = reflash1pt2.Time_sec_;
boostR1pt2        = reflash1pt2.Boost_psi_;
clutchR1pt2       = reflash1pt2.ClutchSw_on_off_;
coolantR1pt2      = reflash1pt2.CoolantTemp_F_;
DAMR1pt2          = reflash1pt2.DynAdvMult_DAM_;
feedbackknR1pt2   = -1 .* reflash1pt2.FeedbackKnock___;
EconomyR1pt2      = reflash1pt2.FuelEconomy_mpg_;
GearR1pt2         = reflash1pt2.GearPosition_gear_;
intakeTR1pt2      = reflash1pt2.IntakeTemp_F_;
MAFR1pt2          = reflash1pt2.MAF_g_s_;
manifoldPresR1pt2 = reflash1pt2.ManAbsPress_psi_;
RPMR1pt2          = reflash1pt2.RPM_RPM_;
throttlePosR1pt2  = reflash1pt2.ThrottlePos___;
speedR1pt2        = reflash1pt2.VehicleSpeed_mph_;
wasteR1pt2        = reflash1pt2.WastegateDuty___;

% Reflash 2
timeR2         = reflash2.Time_sec_;
boostR2        = reflash2.Boost_psi_;
clutchR2       = reflash2.ClutchSw_on_off_;
coolantR2      = reflash2.CoolantTemp_F_;
DAMR2          = reflash2.DynAdvMult_DAM_;
feedbackknR2   = -1 .* reflash2.FeedbackKnock___;
EconomyR2      = reflash2.FuelEconomy_mpg_;
GearR2         = reflash2.GearPosition_gear_;
intakeTR2      = reflash2.IntakeTemp_F_;
MAFR2          = reflash2.MAF_g_s_;
manifoldPresR2 = reflash2.ManAbsPress_psi_;
RPMR2          = reflash2.RPM_RPM_;
throttlePosR2  = reflash2.ThrottlePos___;
speedR2        = reflash2.VehicleSpeed_mph_;
wasteR2        = reflash2.WastegateDuty___;

%% Percents

% Convert Control 1 variables to percent (with time in minutes)
timeCP         = timeC ./ 60;
boostCP        = (boostC ./ max(boostC)) * 100;
%clutchCP       = (clutchC ./ max(clutchC)) * 100;
coolantCP      = (coolantC ./ max(coolantC)) * 100;
DAMCP          = (DAMC ./ max(DAMC)) * 100;
feedbackknCP   = (feedbackknC ./ max(feedbackknC)) * 100;
EconomyCP      = (EconomyC ./ max(EconomyC)) * 100;
GearCP         = (GearC ./ max(GearC)) * 100;
intakeTCP      = (intakeTC ./ max(intakeTC)) * 100;
MAFCP          = (MAFC ./ max(MAFC)) * 100;
manifoldPresCP = (manifoldPresC ./ max(manifoldPresC)) * 100;
RPMCP          = (RPMC ./ max(RPMC)) * 100;
throttlePosCP  = (throttlePosC ./ max(throttlePosC)) * 100;
speedCP        = (speedC ./ max(speedC)) * 100;
wasteCP        = (wasteC ./ max(wasteC)) * 100;

% Convert Control 2 variables to percent
timeC2P         = timeC2 ./ 60;
boostC2P        = (boostC2 ./ max(boostC2)) * 100;
%clutchC2P       = (clutchC2 ./ max(clutchC2)) * 100;
coolantC2P      = (coolantC2 ./ max(coolantC2)) * 100;
DAMC2P          = (DAMC2 ./ max(DAMC2)) * 100;
feedbackknC2P   = (feedbackknC2 ./ max(feedbackknC2)) * 100;
EconomyC2P      = (EconomyC2 ./ max(EconomyC2)) * 100;
GearC2P         = (GearC2 ./ max(GearC2)) * 100;
intakeTC2P      = (intakeTC2 ./ max(intakeTC2)) * 100;
MAFC2P          = (MAFC2 ./ max(MAFC2)) * 100;
manifoldPresC2P = (manifoldPresC2 ./ max(manifoldPresC2)) * 100;
RPMC2P          = (RPMC2 ./ max(RPMC2)) * 100;
throttlePosC2P  = (throttlePosC2 ./ max(throttlePosC2)) * 100;
speedC2P        = (speedC2 ./ max(speedC2)) * 100;
wasteC2P        = (wasteC2 ./ max(wasteC2)) * 100;

% Convert Reflash 1 variables to percent
timeRP         = timeR ./ 60;
boostRP        = (boostR ./ max(boostR)) * 100;
%clutchRP       = (clutchR ./ max(clutchR)) * 100;
coolantRP      = (coolantR ./ max(coolantR)) * 100;
DAMRP          = (DAMR ./ max(DAMR)) * 100;
feedbackknRP   = (feedbackknR ./ max(feedbackknR)) * 100;
EconomyRP      = (EconomyR ./ max(EconomyR)) * 100;
GearRP         = (GearR ./ max(GearR)) * 100;
intakeTRP      = (intakeTR ./ max(intakeTR)) * 100;
MAFRP          = (MAFR ./ max(MAFR)) * 100;
manifoldPresRP = (manifoldPresR ./ max(manifoldPresR)) * 100;
RPMRP          = (RPMR ./ max(RPMR)) * 100;
throttlePosRP  = (throttlePosR ./ max(throttlePosR)) * 100;
speedRP        = (speedR ./ max(speedR)) * 100;
wasteRP        = (wasteR ./ max(wasteR)) * 100;

% Convert Reflash 1 pt 2 variables to percent
timeR1pt2P         = timeR1pt2 ./ 60;
boostR1pt2P        = (boostR1pt2 ./ max(boostR1pt2)) * 100;
coolantR1pt2P      = (coolantR1pt2 ./ max(coolantR1pt2)) * 100;
DAMR1pt2P          = (DAMR1pt2 ./ max(DAMR1pt2)) * 100;
feedbackknR1pt2P   = (feedbackknR1pt2 ./ max(feedbackknR1pt2)) * 100;
EconomyR1pt2P      = (EconomyR1pt2 ./ max(EconomyR1pt2)) * 100;
GearR1pt2P         = (GearR1pt2 ./ max(GearR1pt2)) * 100;
intakeTR1pt2P      = (intakeTR1pt2 ./ max(intakeTR1pt2)) * 100;
MAFR1pt2P          = (MAFR1pt2 ./ max(MAFR1pt2)) * 100;
manifoldPresR1pt2P = (manifoldPresR1pt2 ./ max(manifoldPresR1pt2)) * 100;
RPMR1pt2P          = (RPMR1pt2 ./ max(RPMR1pt2)) * 100;
throttlePosR1pt2P  = (throttlePosR1pt2 ./ max(throttlePosR1pt2)) * 100;
speedR1pt2P        = (speedR1pt2 ./ max(speedR1pt2)) * 100;
wasteR1pt2P        = (wasteR1pt2 ./ max(wasteR1pt2)) * 100;

% Convert Reflash 2 variables to percent
timeR2P         = timeR2 ./ 60;
boostR2P        = (boostR2 ./ max(boostR2)) * 100;
%clutchR2P       = (clutchR2 ./ max(clutchR2)) * 100;
coolantR2P      = (coolantR2 ./ max(coolantR2)) * 100;
DAMR2P          = (DAMR2 ./ max(DAMR2)) * 100;
feedbackknR2P   = (feedbackknR2 ./ max(feedbackknR2)) * 100;
EconomyR2P      = (EconomyR2 ./ max(EconomyR2)) * 100;
GearR2P         = (GearR2 ./ max(GearR2)) * 100;
intakeTR2P      = (intakeTR2 ./ max(intakeTR2)) * 100;
MAFR2P          = (MAFR2 ./ max(MAFR2)) * 100;
manifoldPresR2P = (manifoldPresR2 ./ max(manifoldPresR2)) * 100;
RPMR2P          = (RPMR2 ./ max(RPMR2)) * 100;
throttlePosR2P  = (throttlePosR2 ./ max(throttlePosR2)) * 100;
speedR2P        = (speedR2 ./ max(speedR2)) * 100;
wasteR2P        = (wasteR2 ./ max(wasteR2)) * 100;

%% Plotting For Control 1 & 2

figure('Name','Control Driving Data 1')
% Plot left y-axis variables with unique colors
yyaxis left
h1 = plot(timeCP, boostCP, 'b', 'LineWidth', 1.5, 'Marker', 'none');              % Blue for Boost
hold on; grid on;
h2 = plot(timeCP, RPMCP, 'r', 'LineWidth', 1.5, 'Marker', 'none');                  % Red for RPM
h3 = plot(timeCP, DAMC .* 100, 'g', 'LineWidth', 1.5, 'Marker', 'none');             % Green for DAM
h4 = plot(timeCP, feedbackknCP, 'c', 'LineWidth', 1.5, 'Marker', 'none');            % Cyan for Feedback Knock
h5 = plot(timeCP, GearC .* 10, 'm', 'LineWidth', 1.5, 'Marker', 'none');             % Magenta for Gear
h6 = plot(timeCP, throttlePosCP, 'k', 'LineWidth', 1.5, 'Marker', 'none');           % Black for Throttle Position
h7 = plot(timeCP, wasteCP, 'Color', [1, 0.5, 0], 'LineWidth', 1.5, 'Marker', 'none');% Orange for Waste
% Set left y-axis label
ylabel('Percent of Max');

% Plot right y-axis variable with a unique color
yyaxis right
h8 = plot(timeCP, speedC, 'Color', [0.5, 0, 0.5], 'LineWidth', 1.5, 'Marker', 'none'); % Purple for Speed
ylabel('Speed (mph)');
xlabel('Time (min)')
% Add title and legend
title('Control Driving Data 1: Percent Scaled Data');
legend([h1, h2, h3, h4, h5, h6, h7, h8], ...
       {'Boost', 'RPM', 'DAM', 'Feedback Knock', 'Gear',...
        'Throttle Position', 'Waste', 'Speed'}, ...
       'Location', 'northeastoutside');

%% Control 2

figure('Name','Control Driving Data 2')
% Plot left y-axis variables with unique colors
yyaxis left
h1 = plot(timeC2P, boostC2P, 'b', 'LineWidth', 1.5, 'Marker', 'none');              % Blue for Boost
hold on; grid on;
h2 = plot(timeC2P, RPMC2P, 'r', 'LineWidth', 1.5, 'Marker', 'none');                  % Red for RPM
h3 = plot(timeC2P, DAMC2 .* 100, 'g', 'LineWidth', 1.5, 'Marker', 'none');            % Green for DAM
h4 = plot(timeC2P, feedbackknC2P, 'c', 'LineWidth', 1.5, 'Marker', 'none');           % Cyan for Feedback Knock
h5 = plot(timeC2P, GearC2 .* 10, 'm', 'LineWidth', 1.5, 'Marker', 'none');            % Magenta for Gear
h6 = plot(timeC2P, throttlePosC2P, 'k', 'LineWidth', 1.5, 'Marker', 'none');         % Black for Throttle Position
h7 = plot(timeC2P, wasteC2P, 'Color', [1, 0.5, 0], 'LineWidth', 1.5, 'Marker', 'none');% Orange for Waste
% Set left y-axis label
ylabel('Percent of Max');

% Plot right y-axis variable with a unique color
yyaxis right
h8 = plot(timeC2P, speedC2, 'Color', [0.5, 0, 0.5], 'LineWidth', 1.5, 'Marker', 'none'); % Purple for Speed
ylabel('Speed (mph)');
xlabel('Time (min)')
% Add title and legend
title('Control Driving Data 2: Percent Scaled Data');
legend([h1, h2, h3, h4, h5, h6, h7, h8], ...
       {'Boost', 'RPM', 'DAM', 'Feedback Knock', 'Gear',...
        'Throttle Position', 'Waste', 'Speed'}, ...
       'Location', 'northeastoutside');

%% Reflash 1

figure('Name','Reflash Driving Data 1')
% Plot left y-axis variables with unique colors
yyaxis left
h1 = plot(timeRP, boostRP, 'b', 'LineWidth', 1.5, 'Marker', 'none');              % Blue for Boost
hold on; grid on;
h2 = plot(timeRP, RPMRP, 'r', 'LineWidth', 1.5, 'Marker', 'none');                  % Red for RPM
h3 = plot(timeRP, DAMR .* 100, 'g', 'LineWidth', 1.5, 'Marker', 'none');            % Green for DAM
h4 = plot(timeRP, feedbackknRP, 'c', 'LineWidth', 1.5, 'Marker', 'none');           % Cyan for Feedback Knock
h5 = plot(timeRP, GearR .* 10, 'm', 'LineWidth', 1.5, 'Marker', 'none');            % Magenta for Gear
h6 = plot(timeRP, throttlePosRP, 'k', 'LineWidth', 1.5, 'Marker', 'none');          % Black for Throttle Position
h7 = plot(timeRP, wasteRP, 'Color', [1, 0.5, 0], 'LineWidth', 1.5, 'Marker', 'none');% Orange for Waste
% Set left y-axis label
ylabel('Percent of Max');

% Plot right y-axis variable with a unique color
yyaxis right
h8 = plot(timeRP, speedR, 'Color', [0.5, 0, 0.5], 'LineWidth', 1.5, 'Marker', 'none'); % Purple for Speed
ylabel('Speed (mph)');
xlabel('Time (min)')
% Add title and legend
title('Reflash Driving Data 1: Percent Scaled Data');
legend([h1, h2, h3, h4, h5, h6, h7, h8], ...
       {'Boost', 'RPM', 'DAM', 'Feedback Knock', 'Gear',...
        'Throttle Position', 'Waste', 'Speed'}, ...
       'Location', 'northeastoutside');

%% Reflash 1 pt 2

figure('Name','Reflash Driving Data 1 pt2')
% Plot left y-axis variables with unique colors
yyaxis left
h1 = plot(timeR1pt2P, boostR1pt2P, 'b', 'LineWidth', 1.5, 'Marker', 'none');              % Blue for Boost
hold on; grid on;
h2 = plot(timeR1pt2P, RPMR1pt2P, 'r', 'LineWidth', 1.5, 'Marker', 'none');                  % Red for RPM
h3 = plot(timeR1pt2P, DAMR1pt2 .* 100, 'g', 'LineWidth', 1.5, 'Marker', 'none');             % Green for DAM
h4 = plot(timeR1pt2P, feedbackknR1pt2P, 'c', 'LineWidth', 1.5, 'Marker', 'none');            % Cyan for Feedback Knock
h5 = plot(timeR1pt2P, GearR1pt2 .* 10, 'm', 'LineWidth', 1.5, 'Marker', 'none');             % Magenta for Gear
h6 = plot(timeR1pt2P, throttlePosR1pt2P, 'k', 'LineWidth', 1.5, 'Marker', 'none');           % Black for Throttle Position
h7 = plot(timeR1pt2P, wasteR1pt2P, 'Color', [1, 0.5, 0], 'LineWidth', 1.5, 'Marker', 'none');% Orange for Waste
% Set left y-axis label
ylabel('Percent of Max');

% Plot right y-axis variable with a unique color
yyaxis right
h8 = plot(timeR1pt2P, speedR1pt2, 'Color', [0.5, 0, 0.5], 'LineWidth', 1.5, 'Marker', 'none'); % Purple for Speed
ylabel('Speed (mph)');
xlabel('Time (min)')
% Add title and legend
title('Reflash Driving Data 1 pt2: Percent Scaled Data');
legend([h1, h2, h3, h4, h5, h6, h7, h8], ...
       {'Boost', 'RPM', 'DAM', 'Feedback Knock', 'Gear',...
        'Throttle Position', 'Waste', 'Speed'}, ...
       'Location', 'northeastoutside');

%% Reflash 2

figure('Name','Reflash Driving Data 2')
% Plot left y-axis variables with unique colors
yyaxis left
h1 = plot(timeR2P, boostR2P, 'b', 'LineWidth', 1.5, 'Marker', 'none');              % Blue for Boost
hold on; grid on;
h2 = plot(timeR2P, RPMR2P, 'r', 'LineWidth', 1.5, 'Marker', 'none');                  % Red for RPM
h3 = plot(timeR2P, DAMR2 .* 100, 'g', 'LineWidth', 1.5, 'Marker', 'none');             % Green for DAM
h4 = plot(timeR2P, feedbackknR2P, 'c', 'LineWidth', 1.5, 'Marker', 'none');            % Cyan for Feedback Knock
h5 = plot(timeR2P, GearR2 .* 10, 'm', 'LineWidth', 1.5, 'Marker', 'none');             % Magenta for Gear
h6 = plot(timeR2P, throttlePosR2P, 'k', 'LineWidth', 1.5, 'Marker', 'none');           % Black for Throttle Position
h7 = plot(timeR2P, wasteR2P, 'Color', [1, 0.5, 0], 'LineWidth', 1.5, 'Marker', 'none');% Orange for Waste
% Set left y-axis label
ylabel('Percent of Max');

% Plot right y-axis variable with a unique color
yyaxis right
h8 = plot(timeR2P, speedR2, 'Color', [0.5, 0, 0.5], 'LineWidth', 1.5, 'Marker', 'none'); % Purple for Speed
ylabel('Speed (mph)');
xlabel('Time (min)')
% Add title and legend
title('Reflash Driving Data 2: Percent Scaled Data');
legend([h1, h2, h3, h4, h5, h6, h7, h8], ...
       {'Boost', 'RPM', 'DAM', 'Feedback Knock', 'Gear',...
        'Throttle Position', 'Waste', 'Speed'}, ...
       'Location', 'northeastoutside');

%% Roshan Jaiswal-Ferri
%Section - 01 
%Aero 303 SSWT Calcs: 2/20/25

%% Workspace Prep

format long     %Allows for more accurate decimals
close all;      %Clears all
clear all;      %Clears Workspace
clc;            %Clears Command Window

%% Data Manip

Data = readtable("AERO303_SSWT_DataM.xlsx");

t = Data.DTime_seconds_;
M = Data.Mach;


Data2 = readtable("AERO303_SSWT_Data2M.xlsx");

t2 = Data2.DTime_seconds_;
M2 = Data2.Mach;

figure;
subplot(2,1,1); 
plot(t, M);
grid on;
xlabel('Time (Seconds)');
ylabel('Mach #');
title('Mach vs Time SSWT Run 1');
subplot(2,1,2); 
plot(t2, M2);
grid on;
xlabel('Time (Seconds)');
ylabel('Mach #');
title('Mach vs Time SSWT Run 2');


%%

beta2 = [69,64,57,49,40,38,33,30,28];
x = 1;
g = 1.4;
M = 3.13;

for i = 1:1:length(beta2)
    B = deg2rad(beta2(x));
    %theta(x) = atand(2*cot(B) * (((M^2)*(sin(B)^2))-1)/(((M^2)*(g+cos(2*B)))+2));
    theta(x) = betatotheta(M,B,g);
    x = x + 1;
end

% Define range of Mach numbers
M_values = [2, 3, 4, 5, 6, 8];  % Example Mach numbers
gamma = 1.4;  % Specific heat ratio for air

% Define range of shock wave angles (beta) in degrees
beta = linspace(10, 90, 200);  % Shock angle from 10° to 90°

% Initialize figure
figure;
plot(theta,beta2,'-*','DisplayName','M = 3.13')
hold on
grid on

% Loop through each Mach number
for M = M_values
    theta = zeros(size(beta)); % Initialize theta array

    % Compute theta for each beta
    for i = 1:length(beta)
        beta_rad = deg2rad(beta(i)); % Convert beta to radians
        theta(i) = rad2deg(atan(2 * cot(beta_rad) * ...
                   ((M^2 * sin(beta_rad)^2 - 1) / ...
                   (M^2 * (gamma + cos(2 * beta_rad)) + 2))));
    end

    % Plot theta vs beta for this Mach number
    plot(theta, beta, 'DisplayName', sprintf('M = %.1f', M));
end

% Labels and legend
ylabel('Shock Wave Angle, \beta (degrees)');
xlabel('Flow Deflection Angle, \theta (degrees)');
title('\theta-\beta-M Chart');
legend show;
xlim([0 45]); % Set x-axis limits
ylim([10 90]);  % Set y-axis limits based on max theta
hold off;


%%

% beta = [57,49,40,38,33,30,28];
% x = 1;
% g = 1.4;
% M = 3.13;
% 
% for i = 1:1:7
%     B = deg2rad(beta(x));
%     %theta(x) = atand(2*cot(B) * (((M^2)*(sin(B)^2))-1)/(((M^2)*(g+cos(2*B)))+2));
%     theta(x) = betatotheta(M,B,g);
%     x = x + 1;
% end
% 
% figure
% plot(theta,beta)
% hold on
% grid on
% plot(theta,beta,'*')
% xlabel('Deflection Angle (\theta)')
% ylabel('Shock Angle (\beta)')
% title('\theta - \beta @Mach 3.13 (In Degrees)')

%% Functions

function theta = betatotheta(M, beta, gamma)
    % wants radians
    theta = atan2(2 * cot(beta) * ((M^2 * sin(beta)^2 - 1) / ...
          (M^2 * (gamma + cos(2 * beta)) + 2)), 1);
    theta = rad2deg(theta); % Convert from radians to degrees
end

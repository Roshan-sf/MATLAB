%% Workspace Preparation
format long; % High precision for calculations
close all;  % Close all figures
clear;      % Clear all variables
clc;        % Clear command window

%% Constants and Parameters
% Gravitational parameters (km^3/s^2)
muSun = 1.327e11; % Sun
muEarth = 398600; % Earth
muJupiter = 126686534; % Jupiter

% Radii (km)
Rearth = 6378; % Earth radius
Rjupiter = 71490; % Jupiter radius

% Orbit parameters
RparkEarth = Rearth + 500; % Earth parking orbit
RparkJupiter = Rjupiter + 20000; % Jupiter parking orbit

% Lambert's problem solver tolerance
lambertTol = 1e-7;

% Time tolerances for ODE integration
odeOptions = odeset('RelTol', 1e-8, 'AbsTol', 1e-8);

%% Mission Timeline
% Define departure and arrival dates
departureDate = datetime(2026, 11, 1, 0, 0, 0);
arrivalDates = [
    datetime(2030, 3, 1, 0, 0, 0);
    datetime(2030, 8, 1, 0, 0, 0);
    datetime(2030, 10, 1, 0, 0, 0)
];

% Convert to Julian Dates
JD2000 = juliandate(datetime(2000, 1, 1, 0, 0, 0));
JDdep = juliandate(departureDate);
JDs = juliandate(arrivalDates);

% Compute transfer times in seconds
transferTimes = (JDs - JDdep) * 86400; % Delta t in seconds

%% Orbital Elements
% Departure planet (Earth)
[earthElements] = AERO351planetary_elements2(3, (JDdep - JD2000) / 36525);
[Re, Ve] = coes2rvd_from_elements(earthElements, muSun);

% Arrival planet (Jupiter) for each arrival date
jupiterElements = arrayfun(@(JD) ...
    AERO351planetary_elements2(5, (JD - JD2000) / 36525), ...
    JDs, 'UniformOutput', false);

jupiterRV = cellfun(@(ele) coes2rvd_from_elements(ele, muSun), ...
    jupiterElements, 'UniformOutput', false);

% Calculate R and V vectors for Jupiter at each arrival date
RJ = cell(1, length(JDs));
VJ = cell(1, length(JDs));
for i = 1:length(JDs)
    elements = AERO351planetary_elements2(5, (JDs(i) - JD2000) / 36525);
    [RJ{i}, VJ{i}] = coes2rvd_from_elements(elements, muSun);
end


%% Lambert's Problem Solutions
% Short way (tm = 1) and long way (tm = -1)
transferTypes = [1, -1]; % Short and long transfers
V1 = cell(2, length(transferTimes));
V2 = cell(2, length(transferTimes));

for tmIdx = 1:length(transferTypes)
    tm = transferTypes(tmIdx);
    for i = 1:length(transferTimes)
        [V1{tmIdx, i}, V2{tmIdx, i}] = lambUVBi(Re, RJ{i}, transferTimes(i), tm, muSun, lambertTol);
    end
end

%% Delta-V Calculations
% Compute parking orbit velocity
VparkEarth = sqrt(muEarth / RparkEarth);
VparkJupiter = sqrt(muJupiter / RparkJupiter);

% Initialize delta-v results
deltaV = struct('departure', zeros(2, length(transferTimes)), ...
                'arrival', zeros(2, length(transferTimes)), ...
                'total', zeros(2, length(transferTimes)));

% Calculate delta-v for departure and arrival burns
for tmIdx = 1:length(transferTypes)
    for i = 1:length(transferTimes)
        % Departure burn
        VinfDep = norm(V1{tmIdx, i} - Ve);
        VburnoutEarth = sqrt(VinfDep^2 + (2 * muEarth / RparkEarth));
        deltaV.departure(tmIdx, i) = abs(VburnoutEarth - VparkEarth);

        % Arrival burn
        VinfArr = norm(VJ{i} - V2{tmIdx, i});
        VburnoutJupiter = sqrt(VinfArr^2 + (2 * muJupiter / RparkJupiter));
        deltaV.arrival(tmIdx, i) = abs(VburnoutJupiter - VparkJupiter);

        % Total delta-v
        deltaV.total(tmIdx, i) = deltaV.departure(tmIdx, i) + deltaV.arrival(tmIdx, i);
    end
end

%% Display Results
arrivalLabels = ["March 1, 2030", "August 1, 2030", "October 1, 2030"];
transferLabels = ["Short Way", "Long Way"];

disp('Delta-V Results (km/s):');
for tmIdx = 1:length(transferTypes)
    fprintf('\n%s:\n', transferLabels(tmIdx));
    for i = 1:length(arrivalDates)
        fprintf('%s - Departure: %.3f km/s, Arrival: %.3f km/s, Total: %.3f km/s\n', ...
            arrivalLabels(i), deltaV.departure(tmIdx, i), ...
            deltaV.arrival(tmIdx, i), deltaV.total(tmIdx, i));
    end
end

% Determine the best arrival date
[~, bestIdx] = min(deltaV.total, [], 2);
for tmIdx = 1:length(transferTypes)
    fprintf('\nBest Arrival Date for %s: %s (%.3f km/s)\n', ...
        transferLabels(tmIdx), arrivalLabels(bestIdx(tmIdx)), ...
        deltaV.total(tmIdx, bestIdx(tmIdx)));
end

%% Plotting
% Transfer orbits and planetary orbits
colors = ['r', 'g', 'b'];
for tmIdx = 1:length(transferTypes)
    figure('Name', sprintf('Transfer Orbits - %s', transferLabels(tmIdx)));
    hold on;
    plot3(0, 0, 0, 'y*', 'MarkerSize', 10); % Sun
    plot3(Re(1), Re(2), Re(3), 'co', 'MarkerSize', 15); % Earth

    % Earth and Jupiter orbits
    plot_planet_orbit(Re, Ve, muSun, odeOptions, 'w--');
    cellfun(@(R, V) plot_planet_orbit(R, V, muSun, odeOptions, 'm--'), RJ, VJ);

    % Transfer orbits
    for i = 1:length(transferTimes)
        state0 = [Re; V1{tmIdx, i}];
        tspan = [0, transferTimes(i)];
        [~, T] = ode45(@twobodymotion, tspan, state0, odeOptions, muSun);
        plot3(T(:, 1), T(:, 2), T(:, 3), colors(i), 'LineWidth', 1.5);
        plot3(T(1, 1), T(1, 2), T(1, 3), [colors(i), '*'], 'MarkerSize', 10); % Start
        plot3(T(end, 1), T(end, 2), T(end, 3), [colors(i), 'o'], 'MarkerSize', 10); % End
    end

    xlabel('X (km)');
    ylabel('Y (km)');
    zlabel('Z (km)');
    grid on;
    legend('Sun', 'Earth', 'Earth Orbit', 'Jupiter Orbits', ...
        'March Transfer', 'August Transfer', 'October Transfer');
    title(sprintf('Transfer Orbits - %s', transferLabels(tmIdx)));
    axis equal;
end

%% Helper Functions
function [R, V] = coes2rvd_from_elements(elements, mu)
    [~, ~, R, V] = coes2rvd(elements(1), elements(2), ...
        elements(3), elements(4), elements(5), elements(6), mu);
end

function plot_planet_orbit(R, V, mu, options, style)
    [~, period] = rv2coes(R, V, mu, 0);
    tspan = [0, period];
    [~, trajectory] = ode45(@twobodymotion, tspan, [R; V], options, mu);
    plot3(trajectory(:, 1), trajectory(:, 2), trajectory(:, 3), style, 'LineWidth', 1.5);
end

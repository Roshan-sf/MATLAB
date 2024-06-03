%Roshan Jaiswal-Ferri
%Section - 03
%Aero 300 Lab 8 - Beam Deflections and PDEs: 5/24/24

%% Workspace Prep

format long     %Allows for more accurate decimals
close all;      %Clears all
clear all;      %Clears Workspace
clc;            %Clears Command Window

%% PART 1: 

% Parameters for mode 2
mode = 2; % Second mode
n_pts = 100; % Number of points
plots = 'no'; % Turn off plots

% Run the function and check the frequency
freq = beamvibe(mode, n_pts, plots);
disp(['Frequency: ', num2str(freq), ' Hz']);


function [freq] = beamvibe(mode, n_pts, plots)
    % Constants
    L = pi; % Length of the beam (3.14 meters)
    C = -9.3979e-6; 
    Q = -0.2; % Initial Q value
    dx = L / (n_pts - 1); % Spacing between points
    dt = sqrt(Q * C * dx^4); % Initial time step

    % Initial conditions
    x = linspace(0, L, n_pts); % x coordinates
    y_initial = sin(mode * x); % Initial deflection
    y = zeros(n_pts, 2); % Array to store deflections at two time steps
    y(:,1) = y_initial; % Set initial deflection
    y(:,2) = y_initial + dt * 0; % Initial velocity assumed to be zero

    % Time stepping
    t_total = 0.1; % Total time for simulation
    t_steps = floor(t_total / dt); % Number of time steps

    % Animation setup
    if strcmp(plots, 'yes')
        figure;
        h = plot(x, y(:,1), 'LineWidth', 2);
        ylim([-1, 1]);
        title('Vibration of Pin-Pin Beam');
        xlabel('Position along the beam (m)');
        ylabel('Deflection (m)');
    end

    for j = 2:t_steps
        % Calculate new deflection using central differencing in space
        for i = 3:n_pts-2
            y(i,3) = Q * (y(i+2,2) - 4*y(i+1,2) + 6*y(i,2) - 4*y(i-1,2) + y(i-2,2)) + 2*y(i,2) - y(i,1);
        end
        
        % Boundary conditions
        y(1,3) = 0; % Pin condition at the left end
        y(2,3) = 0; % Pin condition at the left end
        y(end,3) = 0; % Pin condition at the right end
        y(end-1,3) = 0; % Pin condition at the right end

        % Update for next time step
        y(:,1) = y(:,2);
        y(:,2) = y(:,3);

        % Plotting
        if strcmp(plots, 'yes')
            set(h, 'YData', y(:,2));
            drawnow;
        end
    end

    % Calculate frequency of vibration
    period = 2 * dt * t_steps; % Rough estimate of period
    freq = 1 / period; % Frequency in Hz

    % Final plot for verification
    if strcmp(plots, 'yes')
        figure;
        plot(x, y(:,2), 'LineWidth', 2);
        title('Final Deflection of Pin-Pin Beam');
        xlabel('Position along the beam (m)');
        ylabel('Deflection (m)');
    end
end





















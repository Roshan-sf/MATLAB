close all;              % Clears all
clear all;              % Clears workspace
clc;                    % Clears command line

mode = 2;
n_pts = 10;
plots = 'yes';

[freq1,W] = beamvibe(mode, n_pts, plots);
disp(['Frequency for first mode: ', num2str(freq1), ' Hz']);

mode = 2;
n_pts = 100;
plots = 'no';

[freq2] = beamvibe(mode, n_pts, plots);
disp(['Frequency for second mode: ', num2str(freq2), ' Hz']);

function [freq, W] = beamvibe(mode, n_pts, plots)
    % Constants and initial setup
    L = pi; % Length of the beam
    dx = L / (n_pts - 1); % Spatial step size
    C = -9.3979e-6; % Given constant
    Q = -0.2; % Given stability factor
    dt = sqrt(abs(Q * C * dx^4)); % Initial time step
    t_end = 0.1; % Simulation end time
    time_steps = ceil(t_end / dt); % Number of time steps

    % Initialize W
    W = zeros(n_pts, time_steps+1);

    % Initial conditions
    x = linspace(0, L, n_pts);
    W(:,1) = sin(mode * x); % Initial deflection

    % Arrays to store the time history of central point for frequency analysis
    time = linspace(0, t_end, time_steps+1);

    % Time stepping loop
    j = 1;
    for i = 2:n_pts-1
        if i == 2
            W(i,j+1) = 0.5*Q * (W(i+2,j) - 4*W(i+1,j) + 6*W(i,j) - 4*W(i-1,j) - W(i,j)) + W(i,j);
        elseif i == n_pts-1
            W(i,j+1) = 0.5*Q * (-W(i,j) - 4*W(i+1,j) + 6*W(i,j) - 4*W(i-1,j) + W(i-2,j)) + W(i,j);
        else
            W(i,j+1) = 0.5*Q * (-W(i,j) - 4*W(i+1,j) + 6*W(i,j) - 4*W(i-1,j) + W(i-2,j)) + W(i,j);
        end
    end

    for j = 2:time_steps
        for i = 2:n_pts-1
            if i == 2
                W(i,j+1) = Q * (W(i+2,j) - 4*W(i+1,j) + 6*W(i,j) - 4*W(i-1,j) - W(i,j)) + 2*W(i,j) - W(i,j-1);
            elseif i == n_pts-1
                W(i,j+1) = Q * (-W(i,j) - 4*W(i+1,j) + 6*W(i,j) - 4*W(i-1,j) + W(i-2,j)) + 2*W(i,j) - W(i,j-1);
            else
                W(i,j+1) = Q * (W(i+2,j) - 4*W(i+1,j) + 6*W(i,j) - 4*W(i-1,j) + W(i-2,j)) + 2*W(i,j) - W(i,j-1);
            end
        end

        % Apply boundary conditions
        W(1,j+1) = 0; % Pin at the start
        W(n_pts,j+1) = 0; % Pin at the end

        % Plot the beam if required
        if strcmp(plots, 'yes')
            plot(x, W(:,j+1));
            ylim([-1, 1]);
            title(['Time: ', num2str(j * dt), ' seconds']);
            drawnow;
        end
    end

    % Calculating frequency
    cross = 0;
    time1 = 0;
    for i = 4:length(W(2,:))
        if W(2,i)*W(2,i-1) < 0
            if time1 == 0
                time1 = i * dt;
            end
            cross = cross + 1;
            time2 = i * dt;
        end
    end

    freq = (cross) / (2 * abs(time1 - time2));

    % Plot the displacement of the central point over time
    if strcmp(plots, 'yes')
        figure;
        plot(time, W(ceil(n_pts/2), :));
        title('Displacement of the Central Point Over Time');
        xlabel('Time (seconds)');
        ylabel('Displacement');
        grid on;
    end
end

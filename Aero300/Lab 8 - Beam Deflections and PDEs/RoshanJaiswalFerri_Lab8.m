%Roshan Jaiswal-Ferri
%Section - 03
%Aero 300 Lab 8 - Beam Deflections and PDEs: 5/24/24

%% Workspace Prep

format long     %Allows for more accurate decimals
close all;      %Clears all
clear all;      %Clears Workspace
clc;            %Clears Command Window

%% PART 1: 

anim = input('Animation? y/n', 's');

if strcmp(anim,'y')
    mode = 1;
    n_pts = 100;
    plots = 'yes';
    
    [freq1] = beamvibe(mode, n_pts, plots);
    %disp(['Frequency for mode 1: ', num2str(freq1), ' Hz']);
    
    % mode = 2;
    % n_pts = 100;
    % plots = 'no';
    % 
    % [freq2] = beamvibe(mode, n_pts, plots);
    % disp(['Frequency for mode 2: ', num2str(freq2), ' Hz']);
end

freq = zeros(1,20);
mode = linspace(1,20,20);

for i = 1:20
    n_pts = 100;
    plots = 'no';
    freq(1,i) = beamvibe(i, n_pts, plots);
    disp(['Frequency for mode ', num2str(i) ,': ', num2str(freq(1,i)), ' Hz']);
end

figure('name', 'Mode vs Frequency')
semilogy(mode,freq, '*', 'LineWidth', 5);
xlabel('Mode')
ylabel('Frequency (Hz)')
title('Mode vs Frequency')
grid on


function [freq] = beamvibe(mode, n_pts, plots)
    % Constants and initial setup
    L = pi; % Length of the beam
    dx = L / (n_pts - 1); % Spatial step size
    C = -9.3979e-6; % Given constant
    Q = -0.2; % Given stability factor
    dt = sqrt(Q * C * (dx^4)); % Initial time step
    t_end = 0.1; % Simulation end time
    time_steps = ceil(t_end / dt); % Number of time steps

    % Initialize W
    W = zeros(n_pts, time_steps+1);

    % Initial conditions
    x = linspace(0, L, n_pts);
    W(:,1) = sin(mode * x); % Initial deflection

    % Arrays to store the time history of central point for frequency analysis
    time = linspace(0,t_end,time_steps+1);

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

        %animation code, it will look a little fast, to fix set mod(j,1)
        %and n pts to 20 or less (will give better animation but worse
        %data)
        if strcmp(plots, 'yes')
            if mod(j,80) == 0
                plot(x, W(:,j+1));
                ylim([-1, 1]);
                title(['Time: ', num2str(j * dt), ' seconds']);
                drawnow;
            end
        end
    end

    %Calculating frequency
    cross = 0;
    time1 = 0;
    for i = 4:length(W(2,:))
        if W(2,i)*W(2,i-1) < 0
            if time1 == 0
                time1 = i*dt;
            end
            cross = cross + 1;
            time2 = i*dt;
        end
    end

    freq = (cross)/(2 * abs(time1 - time2));

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

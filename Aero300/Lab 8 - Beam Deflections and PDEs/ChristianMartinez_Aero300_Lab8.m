% Christian Martinez
% Aero 300 Lab 8: 5/24/24

close all;                   % clears all 
clear all;                   % clears workspace
clc;                         % clears command line


% Part 1 of AIAA
mode = 2;
n_pts = 20;
plots = 'no';

[freq1] = beamvibe(mode, n_pts, plots);
disp(['Frequency for first mode: ', num2str(freq1), ' Hz']);

mode = 2;
n_pts = 20;
plots = 'yes';

[freq2] = beamvibe(mode, n_pts, plots);
disp(['Frequency for second mode: ', num2str(freq2), ' Hz']);

function [freq] = beamvibe(mode, n_pts, plots)
    % Constants and initial setup
    L = pi; % length of the beam
    dx = L / (n_pts - 1); % spatial step size
    C = -9.3979e-6; % given constant
    Q = -0.2; % given stability factor
    dt = sqrt(Q * C * dx^4); % initial time step
    t_end = 0.1; % simulation end time
    time_steps = ceil(t_end / dt); % number of time steps

    % Initial conditions
    x = linspace(0, L, n_pts);
    W = sin(mode * x); % initial deflection
    W = W; % initial condition for previous time step
    W_next = zeros(1, n_pts); % array for the next time step

    % Arrays to store the time history of central point for frequency analysis
    W_center = zeros(1, time_steps);
    time = zeros(1, time_steps);

    % Time stepping loop
    for t = 1:time_steps
        % Update internal points using finite difference method
        for i = 3:n_pts-2
            W_next(i) = 0.5 * Q * (W(i+2) - 4*W(i+1) + 6*W(i) - 4*W(i-1) + W(i-2)) + W(i);
        end

        % Apply boundary conditions
        W_next(1) = 0; % Pin at the start
        W_next(end) = 0; % Pin at the end

        % Store the central point displacement for frequency analysis
        W_center(t) = W_next(ceil(n_pts/2));
        time(t) = t * dt;

        % Update the displacement arrays
        W_prev = W;
        W = W_next;

        % Plot the beam if required
        if strcmp(plots, 'yes')
            plot(x, W);
            ylim([-1, 1]);
            title(['Time: ', num2str(t * dt), ' seconds']);
            drawnow;
        end
    end

    % Calculate the frequency using FFT
    Y = fft(W_center);
    P2 = abs(Y / time_steps);
    P1 = P2(1:time_steps/2+1);
    P1(2:end-1) = 2*P1(2:end-1);
    f = (1/dt) * (0:(time_steps/2)) / time_steps;

    % Find the peak frequency in the spectrum
    [~, idx] = max(P1);
    freq = f(idx);

    % Plot frequency spectrum if required
    if strcmp(plots, 'yes')
        figure;
        plot(f, P1)
        title('Single-Sided Amplitude Spectrum of W_{center}(t)')
        xlabel('Frequency (Hz)')
        ylabel('|P1(f)|')
    end
end


%count/2)/delta t i = 4:length(2,:))
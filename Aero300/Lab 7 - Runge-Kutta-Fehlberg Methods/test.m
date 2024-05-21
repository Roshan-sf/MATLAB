%Roshan Jaiswal-Ferri
%Section - 03
%Aero 300 Lab 7 - Runge-Kutta-Fehlberg Methods: 5/17/24

%% Workspace Prep

format long     %Allows for more accurate decimals
close all;      %Clears all
clear all;      %Clears Workspace
clc;            %Clears Command Window

%% PART 1: Using Solving y(t) with rfk45
t = [1, pi/3];
y0 = 1;
h = 0.01;
tol = 1e-6;


[t, y] = rkf45(@fx, t, y0, h, tol);

figure
plot(t,y);
grid on
legend('y','yPrime','Location','best')


%% Defining given ODE for part 1

function dydt = fx(t,y)
    dydt = [y(1);(y-t-1).^2+2]; %isolated y double prime, and create function
end

%% Runge-Kutta-Fehlberg Function

function [t, y] = rkf45(fun, tspan, y0, h, rTol)
% RKF45 Solves a first-order ODE using the Runge-Kutta-Fehlberg method
%
% Inputs:
%   fun   - Function handle representing the ODE (dy/dt = fun(t,y))
%   tspan - [t0 tf] Initial and final times
%   y0    - Initial condition
%   h     - Initial step size
%   rTol  - Relative error tolerance
%
% Outputs:
%   t - Vector of time points
%   y - Vector of solution values at corresponding time points

    % Extract initial and final times
    t0 = tspan(1);
    tf = tspan(2);

    % Initialize variables
    t = t0;
    y = y0;

    % Preallocate arrays for efficiency (adjust size as needed)
    maxSteps = 10000;
    t = zeros(maxSteps, 1);
    y = zeros(maxSteps, 1);

    % Set initial values
    t(1) = t0;
    y(1) = y0;

    % Counters for the loop
    i = 1;

    % Loop until reaching the final time
    while t(i) < tf
        % Ensure that the step size doesn't go beyond the final time
        if t(i) + h > tf
            h = tf - t(i);
        end

        % Compute the Runge-Kutta-Fehlberg coefficients
        k1 = h * feval(fun, t(i), y(i));
        k2 = h * feval(fun, t(i) + 1/4*h, y(i) + 1/4*k1);
        k3 = h * feval(fun, t(i) + 3/8*h, y(i) + 3/32*k1 + 9/32*k2);
        k4 = h * feval(fun, t(i) + 12/13*h, y(i) + 1932/2197*k1 - 7200/2197*k2 + 7296/2197*k3);
        k5 = h * feval(fun, t(i) + h, y(i) + 439/216*k1 - 8*k2 + 3680/513*k3 - 845/4104*k4);
        k6 = h * feval(fun, t(i) + 1/2*h, y(i) - 8/27*k1 + 2*k2 - 3544/2565*k3 + 1859/4104*k4 - 11/40*k5);

        % Compute the 4th and 5th order estimates
        y4 = y(i) + 25/216*k1 + 1408/2565*k3 + 2197/4104*k4 - 1/5*k5;
        y5 = y(i) + 16/135*k1 + 6656/12825*k3 + 28561/56430*k4 - 9/50*k5 + 2/55*k6;

        % Estimate the error
        err = abs(y5 - y4);

        % Check if the error is within the tolerance
        if err <= rTol * abs(y5)
            % Accept the step
            i = i + 1;
            t(i) = t(i-1) + h;
            y(i) = y5;
        end

        % Adjust the step size
        h = h * min(max(0.84 * (rTol * abs(y5) / err)^0.25, 0.1), 4.0);
    end

    % Trim arrays to the actual size
    t = t(1:i);
    y = y(1:i);
end
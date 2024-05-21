% Christian Martinez
% Aero 300 Lab 7: 5/17/24

close all;                   % clears all 
clear all;                   % clears workspace
clc;                         % clears command line

%% Part 1

% Define the differential equation
fun = @(t, y) (y - t - 1)^2 + 2;

% Define parameters
tspan = [0 pi/3];
y0 = 1;
h = 0.01;
rTol = 1e-6;

% Solve using the custom RKF45 function
[t_rkf, y_rkf] = rkf45(fun, tspan, y0, h, rTol);

% Solve using MATLAB's ode45 function
options = odeset('RelTol', rTol, 'AbsTol', 1e-6);
[t_ode, y_ode] = ode45(fun, tspan, y0, options);

% Exact solution
y_exact = @(t) tan(t) + t + 1;

% Initialize error arrays
error_rkf = zeros(size(t_rkf));
error_ode = zeros(size(t_ode));

% Calculate the error for RKF45 using a for loop
for i = 1:length(t_rkf)
    % Calculate the exact solution at the current t_rkf time point
    y_exact_val = y_exact(t_rkf(i));
    
    % Calculate the error for RKF45
    error_rkf(i) = abs(y_rkf(i) - y_exact_val);
end

% Calculate the error for ode45 using a for loop
for i = 1:length(t_ode)
    % Calculate the exact solution at the current t_ode time point
    y_exact_val = y_exact(t_ode(i));
    
    % Calculate the error for ode45
    error_ode(i) = abs(y_ode(i) - y_exact_val);
end

% Plot the results
figure;
plot(t_rkf, y_rkf, 'b');
hold on;
plot(t_ode, y_ode, 'r--');
plot(t_rkf, y_exact(t_rkf), 'k:');
grid on;
xlabel('Time t');
ylabel('Solution y');
legend('RKF45','ode45','Exact Solution')
title('Solution of the ODE using RKF45 and ode45');

% % Plot the errors
% figure;
% plot(t_rkf, error_rkf, 'b', 'DisplayName', 'Error in RKF45');
% hold on;
% plot(t_ode, error_ode, 'r--', 'DisplayName', 'Error in ode45');
% grid on;
% xlabel('Time t');
% ylabel('Error');
% legend;
% title('Error Comparison of RKF45 and ode45');

% Plot Error
figure;
subplot(2,1,1)
plot(t_rkf, error_rkf, 'b');
grid on;
xlabel('Time t');
ylabel('Error');
legend('Error in RKF45')
title('Error of RKF45');

% Plot Error
subplot(2,1,2)
plot(t_ode, error_ode, 'r');
grid on;
xlabel('Time t');
ylabel('Error');
legend('Error in ode45')
title('Error of ode45');

%% Part 2

% Initial conditions
tspan_lorenz = [0, 50];
y0_lorenz = [0; 1; 20];
h = 1;
rTol = 1e-6;

% tspan_lorenz = [0, 50];
% y0_lorenz = [1; 1; 1]; % Initial conditions

% Solve using the custom RKF45 function
[t_rkf_lorenz, y_rkf_lorenz] = rkf45(@lorenz, tspan_lorenz, y0_lorenz, h, rTol);

% Solve using MATLAB's ode45 function
% [t_ode_lorenz, y_ode_lorenz] = ode45(@(t, y) lorenz_equations(t, y, sig, rho, beta), tspan_lorenz, y0_lorenz, options);

% Plot the results in 3D
figure;
plot3(y_rkf_lorenz(:,1), y_rkf_lorenz(:,2), y_rkf_lorenz(:,3), 'b');
hold on;
% plot3(y_ode_lorenz(:,1), y_ode_lorenz(:,2), y_ode_lorenz(:,3), 'r--');
grid on;
xlabel('x');
ylabel('y');
zlabel('z');
legend('RKF45', 'ode45');
title('Solution of the Lorenz Equations');

% Runge-Kutta Function
function [t, y] = rkf45(fun, tspan, y0, h, rTol)
    % RKF45 Solves ODE using Runge-Kutta-Fehlberg method
    %   [t, y] = rkf45(fun, tspan, y0, h, rTol) uses the Runge-Kutta-Fehlberg
    %   method to solve the ODE defined by 'fun' over the interval 'tspan'
    %   with initial condition 'y0', initial step size 'h', and relative
    %   tolerance 'rTol'.
    %
    % Input:
    %   fun - function handle defining the ODE (dy/dt = f(t, y))
    %   tspan - [t0 tf] time interval
    %   y0 - initial condition
    %   h - initial step size
    %   rTol - relative tolerance for error control
    %
    % Output:
    %   t - vector of time points
    %   y - solution vector

    % Initialize
    t0 = tspan(1);
    tf = tspan(2);
    t = t0;
    y = y0;
    
    % Coefficients for the Runge-Kutta-Fehlberg method
    a2 = 1/4;
    a3 = 3/8;
    a4 = 12/13;
    a5 = 1;
    a6 = 1/2;
    
    b21 = 1/4;
    b31 = 3/32;    b32 = 9/32;
    b41 = 1932/2197; b42 = -7200/2197; b43 = 7296/2197;
    b51 = 439/216;  b52 = -8;      b53 = 3680/513; b54 = -845/4104;
    b61 = -8/27;    b62 = 2;      b63 = -3544/2565; b64 = 1859/4104; b65 = -11/40;
    
    r1 = 25/216;   r3 = 1408/2565; r4 = 2197/4104; r5 = -1/5;
    r1z = 16/135;  r3z = 6656/12825; r4z = 28561/56430; r5z = -9/50; r6z = 2/55;
    
    while t(end) < tf
        % Current values
        ti = t(end);
        yi = y(end);

        % Calculate the slopes
        s1 = fun(ti, yi);
        s2 = fun(ti + a2*h, yi + b21*h*s1);
        s3 = fun(ti + a3*h, yi + b31*h*s1 + b32*h*s2);
        s4 = fun(ti + a4*h, yi + b41*h*s1 + b42*h*s2 + b43*h*s3);
        s5 = fun(ti + a5*h, yi + b51*h*s1 + b52*h*s2 + b53*h*s3 + b54*h*s4);
        s6 = fun(ti + a6*h, yi + b61*h*s1 + b62*h*s2 + b63*h*s3 + b64*h*s4 + b65*h*s5);

        % Calculate the 4th and 5th order estimates
        yi_next = yi + h*(r1*s1 + r3*s3 + r4*s4 + r5*s5);
        zi_next = yi + h*(r1z*s1 + r3z*s3 + r4z*s4 + r5z*s5 + r6z*s6);

        % Error estimation
        e = abs(zi_next - yi_next);

        % Relative error test
        if e / abs(yi_next) < rTol
            % Accept step
            t = [t; ti + h];
            y = [y; yi_next];
        end

        % Adjust step size
        h = 0.8 * (rTol * abs(yi_next) / e)^(1/5) * h;

        % Ensure h doesn't exceed the interval
        if t(end) + h > tf
            h = tf - t(end);
        end
    end
end

% Lorenz Function
function yp = lorenz(t, y)
    % Parameters
    sig = 10;
    rho = 28;
    beta = 8 / 3;

    % Initialize yp as a column vector
    yp = zeros(3, 1);

    % Lorenz equations
    yp(1) = sig * (y(2) - y(1));
    yp(2) = y(1) * (rho - y(3)) - y(2);
    yp(3) = y(1) * y(2) - beta * y(3);
end
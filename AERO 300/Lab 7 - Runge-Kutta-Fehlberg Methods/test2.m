% Your Name
% Runge-Kutta-Fehlberg Methods for Solving Differential Equations

%% Workspace Preparation

format long;
close all;
clear;
clc;

%% PART 1: Solving ODE with RKF45 and ode45

% Define the function for the given ODE
fun = @(t, y) (y - t - 1)^2 + 2;

% Define parameters
tspan = [0, pi/3];
y0 = 1;
h = 0.01;
rTol = 1e-6;
h_min = 1e-10; % Minimum step size to prevent infinite loop

% Solve using custom RKF45 method
[t_rkf, y_rkf] = rkf45(fun, tspan, y0, h, rTol, h_min);

% Solve using MATLAB's ode45
opts = odeset('RelTol', rTol, 'AbsTol', rTol);
[t_ode45, y_ode45] = ode45(fun, tspan, y0, opts);

% Exact solution
exact_solution = @(t) tan(t) + t + 1;
y_exact = arrayfun(exact_solution, t_rkf);

% Plot the results
figure('Name', 'ODE Solution Comparison');
plot(t_rkf, y_rkf, 'b', 'LineWidth', 1.5);
hold on;
plot(t_ode45, y_ode45, 'r--', 'LineWidth', 1.5);
plot(t_rkf, y_exact, 'k:', 'LineWidth', 1.5);
legend('RKF45', 'ode45', 'Exact Solution', 'Location', 'Best');
xlabel('Time t');
ylabel('y(t)');
title('Solution of ODE using RKF45, ode45, and Exact Solution');
grid on;

% Plot the errors
figure('Name', 'ODE Solution Errors');
error_rkf = abs(y_rkf - y_exact');
error_ode45 = abs(interp1(t_ode45, y_ode45, t_rkf) - y_exact');
plot(t_rkf, error_rkf, 'b', 'LineWidth', 1.5);
hold on;
plot(t_rkf, error_ode45, 'r--', 'LineWidth', 1.5);
legend('Error RKF45', 'Error ode45', 'Location', 'Best');
xlabel('Time t');
ylabel('Error');
title('Errors of RKF45 and ode45 Compared to Exact Solution');
grid on;

%% PART 2: Solving Lorenz System with RKF45 and ode45

% Define the Lorenz system
lorenz = @(t, y) [10*(y(2) - y(1)); y(1)*(28 - y(3)) - y(2); y(1)*y(2) - 8/3*y(3)];

% Define initial conditions and time span
y0_lorenz = [1; 1; 1];
tspan_lorenz = [0, 50];

% Solve using custom RKF45 method
[t_rkf_lorenz, y_rkf_lorenz] = rkf45(@(t, y) lorenz(t, y), tspan_lorenz, y0_lorenz, h, rTol, h_min);

% Solve using MATLAB's ode45
[t_ode45_lorenz, y_ode45_lorenz] = ode45(@(t, y) lorenz(t, y), tspan_lorenz, y0_lorenz, opts);

% Plot the results in 3D
figure('Name', 'Lorenz System Solution');
plot3(y_rkf_lorenz(:,1), y_rkf_lorenz(:,2), y_rkf_lorenz(:,3), 'b', 'LineWidth', 1.5);
hold on;
plot3(y_ode45_lorenz(:,1), y_ode45_lorenz(:,2), y_ode45_lorenz(:,3), 'r--', 'LineWidth', 1.5);
legend('RKF45', 'ode45', 'Location', 'Best');
xlabel('x');
ylabel('y');
zlabel('z');
title('Lorenz System Solution using RKF45 and ode45');
grid on;

%% Function: Runge-Kutta-Fehlberg Method

function [t, y] = rkf45(fun, tspan, y0, h, rTol, h_min)
    t = tspan(1);
    y = y0(:)';
    tn = tspan(1);
    yn = y0(:)';
    
    % Maximum iterations to prevent infinite loop
    max_iter = 10000;
    iter = 0;
    
    while tn < tspan(2)
        iter = iter + 1;
        if iter > max_iter
            warning('Maximum number of iterations reached. Stopping integration.');
            break;
        end
        
        % Ensure the step does not exceed the end of the interval
        h = min(h, tspan(2) - tn);
        
        % Runge-Kutta-Fehlberg steps
        s1 = fun(tn, yn);
        s2 = fun(tn + 1/4*h, yn + 1/4*h*s1);
        s3 = fun(tn + 3/8*h, yn + 3/32*h*s1 + 9/32*h*s2);
        s4 = fun(tn + 12/13*h, yn + 1932/2197*h*s1 - 7200/2197*h*s2 + 7296/2197*h*s3);
        s5 = fun(tn + h, yn + 439/216*h*s1 - 8*h*s2 + 3680/513*h*s3 - 845/4104*h*s4);
        s6 = fun(tn + 1/2*h, yn - 8/27*h*s1 + 2*h*s2 - 3544/2565*h*s3 + 1859/4104*h*s4 - 11/40*h*s5);
        
        w_next = yn + h*(25/216*s1 + 1408/2565*s3 + 2197/4104*s4 - 1/5*s5);
        z_next = yn + h*(16/135*s1 + 6656/12825*s3 + 28561/56430*s4 - 9/50*s5 + 2/55*s6);
        
        e = abs(z_next - w_next);
        rel_error = e ./ (abs(w_next) + eps);
        
        if all(rel_error < rTol)
            % Accept the step
            tn = tn + h;
            yn = w_next;
            t = [t; tn];
            y = [y; yn];
        end
        
        % Adjust step size
        if any(rel_error == 0)
            h = h * 2; % Increase step size if error is zero
        else
            h = 0.8 * h * min(max((rTol ./ rel_error).^(1/5), 0.1), 5);
        end
        
        % Ensure step size is not too small
        if h < h_min
            warning('Step size too small. Stopping integration.');
            break;
        end
    end
end

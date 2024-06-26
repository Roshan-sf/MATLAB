% Isai Villanueva
% AERO 300
% May 27th, 2024
% Lab 7 - Runge-Kutta-Fehlberg Methods

close;
clear;
clc;

%% Runge-Kutta-Fehlberg (1st Order Diff Equation)

% Define ode45 Variables (as listed from lab)
tspan = [0 pi/3];
y0 = 1;

% Ode function
[to,yo] = ode45(@func, tspan, y0);

% RKF45 Variables (as listed from lab)
h = 0.01;
rTol = 10^-6;

% Rkf45 function
[tr,yr] = rkf45(@fun, tspan, y0, h, rTol);

% Plot for ode45, rkf45, and y(t)
j = linspace(0, pi/3, 100); % Represents t in y(t) = tan(t) + t + 1
g = tan(j) + j + 1;

figure()
plot(to,yo,'o') % ode45 command
grid on
hold on
plot(tr,yr,'xr') %rkf45 command
plot(j, g) % actual function
xlabel('t (Input)')
ylabel('y (Output)')
title('ode45 vs. rkf45 vs. y(t)')
legend('ode45 command','rkf45 command','original y(t)','Location','best')



%% Lorenz Equations

time = [0 50];
y0 = [1; 1; 1];
h = 0.01;
rTol = 1*10^-6;

% ode
[to2, yo2] = ode45(@lorenzf, time, y0);

% rkf
[tr2, yr2] = rkf45(@lorenzf, time, y0, h, rTol);

figure()
plot3(yo2(:,1),yo2(:,2),yo2(:,3)) % ode45 command
grid on
hold on
plot3(yr2(:,1),yr2(:,2),yr2(:,3))% rkf45 command
xlabel("x")
ylabel("y")
zlabel("z")
title("Lorenz Function: ODE45 vs. RKF45")
legend("ode45","rkf45")



%% Functions

% rkf45 function
function [tr,yr] = rkf45(fun, tspan, y0, h, rTol)
    % Redefine tspan, y0
    t = tspan(1);
    y = y0;
    tr = t;
    yr = y';
    while t < tspan(2)
        % Calculate the s values
        s1 = fun(t, y);
        s2 = fun(t + 1/4*h, y + 1/4*h*s1);
        s3 = fun(t + 3/8*h, y + 3/32*h*s1 + 9/32*h*s2);
        s4 = fun(t + 12/13*h, y + 1932/2197*h*s1 - 7200/2197*h*s2 + 7296/2197*h*s3);
        s5 = fun(t + h, y + 439/216*h*s1 - 8*h*s2 + 3680/513*h*s3 - 845/4104*h*s4);
        s6 = fun(t + 1/2*h, y - 8/27*h*s1 + 2*h*s2 - 3544/2565*h*s3 + 1859/4104*h*s4 - 11/40*h*s5);
    
        % Calculate the fourth and fifth y values
        y4 = y + 25/216*h*s1 + 1408/2565*h*s3 + 2197/4104*h*s4 - 1/5*h*s5;
        y5 = y + 16/135*h*s1 + 6656/12825*h*s3 + 28561/56430*h*s4 - 9/50*h*s5 + 2/55*h*s6;
        
        % Calculate the error from y4 and y5
        err = norm(y5 - y4);

        % Consider tolerance
        if err <= rTol
            t = t + h;
            y = y5;

            % Store results
            tr(end+1) = t;
            yr(end+1,:) = y';
        else
            % Reduce step size
            h = 0.8*((rTol/err)^0.2)*h; 
            t = t + h;
        end
    end
end

% ode45 Function
function dydt = func(to,yo)
    dydt = (yo - to - 1)^2 + 2; 
end

% ode for rkf45 function
function dydt = fun(tr,yr)
    dydt = (yr - tr - 1)^2 + 2; 
end

% Lorenz Function
function dydt = lorenzf(t, y)
    % Parameters
    sigma = 10;
    rho = 28;
    beta = 8/3;

    % ODEs
    dydt = [
        sigma * (y(2) - y(1));
        y(1) * (rho - y(3)) - y(2);
        y(1) * y(2) - beta * y(3)
        ];
end
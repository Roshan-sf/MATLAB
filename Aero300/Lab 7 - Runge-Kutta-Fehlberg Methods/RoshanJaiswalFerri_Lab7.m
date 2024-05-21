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
h_min = 1e-10; % Minimum step size to prevent infinite loop

%Given ODE
fun = @(t, y) (y - t - 1)^2 + 2;

[t, y] = rkf45(fun, t, y0, h, tol);
options = odeset('RelTol', tol, 'AbsTol', tol);
[t2, y2] = ode45(fun, [0 pi/3], 1, options);

%% Runge-Kutta-Fehlberg Function

function [t, Wn] = rkf45(fun, tspan, y0, h, rTol)
    f = fun;
    y = y0;
    t = tspan(1);
    Wn = y0;
    Zn = y0;
    maxk = 10000;
    n = 1;
    k = 1; %iteration reset

    while t(k) < tspan(1,2)
        if k == maxk
            warning('More than 10,000 Iterations');
            break
        end

        s1 = f(t(k), Wn(k,:));
        s2 = f(t(k) + (1/4)*h, Wn(k,:) + (1/4)*h*s1);
        s3 = f(t(k) + (3/8)*h, Wn(k,:) + (3/32)*h*s1 + (9/32)*h*s2);
        s4 = f(t(k) + (12/13)*h, Wn(k,:) + (1932/2197)*h*s1 - (7200/2197)*h*s2 + (7296/2197)*h*s3);
        s5 = f(t(k) + h, Wn(k,:) + (439/216)*h*s1 - 8*h*s2 + (3680/513)*h*s3 - (845/4104)*h*s4);
        s6 = f(t(k) + (1/2)*h, Wn(k,:) - (8/27)*h*s1 + 2*h*s2 - (3544/2565)*h*s3 + (1859/4104)*h*s4 - (11/40)*h*s5);

        Wn(k+1,:) = Wn(k,:) + h*((25/216)*s1 + (1408/2565)*s3 + (2197/4104)*s4 - (1/5)*s5);
        Zn(k+1,:) = Zn(k,:) + h*((16/135)*s1 + (6656/12825)*s3 + (28561/56430)*s4 - (9/50)*s5 + (2/55)*s6);

        err = norm(Zn(k+1,:)-Wn(k+1,:));

        if (err/norm(Wn(k+1,:))) < rTol
            t(k+1) = t(k) + h;
            Wn(k+1,:) = Zn(k+1,:);
            h = 0.8*(((rTol*norm(Wn(k+1,:)))/(err))^(1/5))*h;
            k = k + 1;
        else
            h = 0.8*(((rTol*norm(Wn(k+1,:)))/(err))^(1/5))*h;
            n = n + 1;
            % if n == 10000
            %     warning('h too small/too many loops')
            %     break
            % end
        end


    end
end







%% Roshan Jaiswal-Ferri
%Section - 02
%Aero 331 HW 3: 3/14/25

%% Workspace Prep

warning off
format long     %Allows for more accurate decimals
close all;      %Clears all
clear all;      %Clears Workspace
clc;            %Clears Command Window

%% Problem 2

A = [1e-4, 2e-4, 3e-4]; % Areas in m^2
b = 0.2; % m
h = 0.4; % m
d = 0.2; % m
G = 26e9; % Pa
t = 0.001; % m
Vz = 10000; % N

a = b * h; % m^2
Yc = (-20 + 3*20) / 6; % cm
Iyy = sum([1, 2, 3] .* (20^2) * 2); % m^4

syms q [1 7]

eqs = [
    q(2) == q(1) - 2083.33;
    q(3) == q(1);
    q(5) == q(4) + 6250;
    q(6) == q(4);
    q(7) == q(6) - q(1) - 4166;
    q(4) == 4166 + q(7) + q(3);
    Vz * (d + b) == (q(1)*b)*h + (q(5)*h)*2*b + (q(6)*b)*h - (q(7)*h)*b;
    (q(1)*b + q(2)*h + q(3)*b - q(7)*h) == (q(6)*b + q(7)*h + q(4)*b + q(5)*h)
];

%solve for q values
soln = solve(eqs, q);
q_vals = double(struct2array(soln));

%angular deformations
a1 = (q_vals(1)*b + q_vals(2)*h + q_vals(3)*b - q_vals(7)*h) / (2*G*a*t); % rads
a2 = (q_vals(6)*b + q_vals(7)*h + q_vals(4)*b - q_vals(5)*h) / (2*G*a*t); % rads

disp('q values:');
disp(num2str(q_vals'));
disp(['a1: ', num2str(a1), ' rads']);
disp(['a2: ', num2str(a2), ' rads']);


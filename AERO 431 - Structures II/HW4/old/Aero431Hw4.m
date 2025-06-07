clc
close all
clear all
%% Problem 2
% Part 1
syms a % crack length [m]

w = 0.5; % specimen width [m]
sigma = 50e6; % applied stress [Pa]

% Define alpha = a/W
alpha = a / w;

% Define the beta function (geometry factor)
beta = (1.122 - 1.122*alpha - 0.820*alpha^2 + 3.768*alpha^3 - 3.040*alpha^4) / sqrt(1 - 2*alpha);

% Stress intensity factor K_I
K1 = beta * sigma * sqrt(pi * a); % [Pa·sqrt(m)]

K1c = 24e6; % MPa*sqrt(m) 

eqn = K1 == K1c;

a = min(real(double(solve(eqn, a))));
disp("Critical Crack Length: " + a +" m")

% Part 2
m = 3.59;
C = 3.15e-11;
sig_max = 50e6;       % Maximum stress [Pa]

a0 = a/2;             % Initial crack length [m]
ac = a0*2;

dsig = 50; % MPa (keep these units for emirical relation)
dN = 100;

Ncycles = 0; % Initialize
Ncycles2 = 0;
a = a0;
while a < ac
    alpha1 = a / w;
    beta1 = (1.122 - 1.122*alpha1- 0.820*alpha1^2 + 3.768*alpha1^3 - 3.040*alpha1^4) / sqrt(1 - 2*alpha1);
    dK = beta1*dsig*sqrt(pi*a);

    a = a + C*dN*dK^m;
    Ncycles = Ncycles + dN;
end

disp("Number of Cycles(dN = 100): " + Ncycles)

dN = 10;
a = a0;
while a < ac
    alpha1 = a / w;
    beta1 = (1.122 - 1.122*alpha1- 0.820*alpha1^2 + 3.768*alpha1^3 - 3.040*alpha1^4) / sqrt(1 - 2*alpha1);
    dK = beta1*dsig*sqrt(pi*a);
    a = a + C*dN*dK^m;
    Ncycles2 = Ncycles2 + dN;
end

disp("Number of Cycles(dN = 10): " + Ncycles2)
%% Comment:
% The difference between the dN = 100 and dN = 10 cycles is not very large.
% (only about 80 cycles, which is negligible)
%% Roshan Jaiswal-Ferri
%Aero 452 Homework 1: 9/24/25

%% Workspace Prep

format long     %Allows for more accurate decimals
close all;      %Clears all
clear all;      %Clears Workspace
clc;            %Clears Command Window

%% PART 1: 

% Known parameters
Ae_At = 5.44;   % Area ratio (Ae/At)
gamma = 1.4;    % Specific heat ratio (can set to 1.2 or 1.4)

% Define symbolic variable
syms Me positive

% Isentropic area-Mach number relation
eqn = Ae_At == (1/Me) * ((2/(gamma+1)) * (1 + ((gamma-1)/2)*Me^2))^((gamma+1)/(2*(gamma-1)));

% Solve numerically for Me (supersonic solution)
Me_sol = vpasolve(eqn, Me, [1 10]);

% Display result
fprintf('Exit Mach number (gamma = %.2f) = %.3f\n', gamma, double(Me_sol));

%%

% Given parameters
g   = 1.3;          % gamma
eps = 6;            % expansion ratio A_e/A_t
pc  = 20;           % chamber pressure [atm]
paSL  = 1.0;        % sea-level ambient pressure [atm]
pa25  = 0.0248;     % ambient pressure at ~25 km [atm] (standard atmosphere)

% Solve area–Mach relation symbolically for the supersonic root
syms M positive
area_eqn = eps == (1/M) * ((2/(g+1))*(1 + (g-1)/2*M^2))^((g+1)/(2*(g-1)));
Me = vpasolve(area_eqn, M, 3);   % initial guess near 3 for supersonic solution

% Exit pressure ratio
pe_pc = (1 + (g-1)/2 * Me^2)^(-g/(g-1));

% Ideal thrust-coefficient core term (pressure-independent part)
CF_core = sqrt((2*g^2)/(g-1) * (2/(g+1))^((g+1)/(g-1)) * (1 - (pe_pc)^((g-1)/g)));

% Add pressure-thrust contribution: (p_e - p_a)/p_c * eps
CF_SL  = CF_core + (pe_pc - paSL/pc)*eps;
CF_25  = CF_core + (pe_pc - pa25/pc)*eps;

% Percent change in thrust (same as change in C_F since F = C_F * p_c * A_t)
pct = double((CF_25 - CF_SL)/CF_SL * 100);

% Display
fprintf('M_e          = %.6f\n', double(Me));
fprintf('p_e/p_c      = %.6f\n', double(pe_pc));
fprintf('C_F (SL)     = %.6f\n', double(CF_SL));
fprintf('C_F (25 km)  = %.6f\n', double(CF_25));
fprintf('Thrust change from SL to 25 km = %.3f %%\n', pct);
    

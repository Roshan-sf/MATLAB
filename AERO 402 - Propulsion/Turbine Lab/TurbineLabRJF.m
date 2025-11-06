%% Roshan Jaiswal-Ferri
%Aero 402 Lab 2: 10/28/25

%% Workspace Prep

warning off
format long     %Allows for more accurate decimals
close all;      %Clears all
clear all;      %Clears Workspace
clc;            %Clears Command Window

%% Plots 1-3

% ============================================================
% PROPS LAB 2 DATA
% ============================================================
pctRPM    = [67 73 80 85 90 95 100];
T0i_F     = [62 63 61 61 62 61 60];
Pi_psig   = [-0.435 -0.545 -0.710 -0.832 -1.009 -1.190 -1.342];
P0i_psig  = [-0.034 -0.044 -0.046 -0.059 -0.065 -0.076 -0.079];
P03_psig  = [17.97 22.22 27.63 31.46 36.61 42.20 47.02];
P07_psig  = [0.27, 0.359, 0.484, 0.426, 0.397, 0.349, 0.344];
P7_psig   = [-0.027, -0.035, -0.091, -0.169, -0.269, -0.368, -0.437];

% ============================================================
% AMBIENT CONDITIONS (from header)
% ============================================================
Ta_F            = 62.2;
Pa_inHg         = 29.92;
psi_per_inHg    = 14.696 / 29.92;
Pamb_psia       = Pa_inHg * psi_per_inHg;
psf_per_psi     = 144;
R_air           = 1716;
gamma           = 1.4;
T_std_R         = 519.67;
P_std_psia      = 14.696;
R_cold          = 1716; %slugs
R_hot           = 1544; 

% ============================================================
% BELLMOUTH INLET GEOMETRY
% ============================================================
ID_in   = 5.5;
r_ft    = (ID_in / 2) / 12;
A_in    = pi * r_ft^2;

% ============================================================
% PRESSURES AND TEMPERATURES (absolute)
% ============================================================
T0i_R    = T0i_F + 459.67;               % inlet stagnation temp, R
Pi_psia  = Pi_psig  + Pamb_psia;         % static absolute, psia
P0i_psia = P0i_psig + Pamb_psia;         % total absolute, psia
P03_psia = P03_psig + Pamb_psia;         % total absolute, psia
P07_psia = P07_psig + Pamb_psia;
P7_psia = P7_psig + Pamb_psia;

Pamb_psf = Pamb_psia * psf_per_psi;
P7_psf = P7_psia * psf_per_psi;
P07_psf = P07_psia * psf_per_psi;
Pi_psf   = Pi_psia  * psf_per_psi;       % psf
P0i_psf  = P0i_psia * psf_per_psi;       % psf


% ============================================================
% INLET FLOW CALCULATIONS
% ============================================================
PR0_static = P0i_psf ./ Pi_psf;
M_in       = sqrt( (2 / (gamma - 1)) .* (PR0_static.^((gamma - 1) / gamma) - 1) );

% Static temperature at inlet
Ti_R = T0i_R ./ (1 + (gamma - 1) / 2 .* M_in.^2);

% Speed of sound and velocity
a_in = sqrt(gamma * R_air .* Ti_R);      % ft/s
V_in = M_in .* a_in;                     % ft/s

% Density from static state
rho_in = Pi_psf ./ (R_air .* Ti_R);

% Mass flow
mdot   = rho_in .* V_in .* A_in;
m_corr = mdot .* sqrt(Ti_R ./ T_std_R) ./ (Pi_psia ./ P_std_psia);
PR_comp = P03_psia ./ P0i_psia;

% ============================================================
% COMPRESSOR MAP PLOT
% ============================================================
figure; hold on; box on; grid on;
plot(m_corr, PR_comp, 'o-', 'LineWidth', 1.5, 'MarkerSize', 7);
xlabel('Corrected Air Flow (slug/s)');
ylabel('Compressor Pressure Ratio P_{03}/P_{01}');
title('Compressor PR vs Corrected Air Flow (bellmouth data only)');
text(m_corr, PR_comp, compose(' %d%%', pctRPM), 'VerticalAlignment','bottom');

% ============================================================
% DATA TABLE
% ============================================================
T = table(pctRPM.', M_in.', mdot.', m_corr.', PR_comp.', ...
    'VariableNames', {'PctRPM','M_in','mdot_slugps','mdot_corr_slugps','PRc'});
disp(T);

% ============================================================
% COMPRESSOR EFFICIENCY
% ============================================================
gamma_g = 1.33; %hot exhaust
R_g     = 1716;

% Ambient for thrust corrections
Ta_R        = Ta_F + 459.67;
T_std_R     = 519.67;
P_std_psia  = 14.696;
Pa_psia     = Pamb_psia;
Pa_psf      = Pa_psia * psf_per_psi;

% Jet-A fuel density
rho_f_lbgal   = 6.7;
g_c           = 32.174;
slug_per_lbm  = 1 / g_c;

% Temperatures
T03_F = [197 247 282 314 342 366 396];
T03_R = T03_F + 459.67;

P03_psf = P03_psia * psf_per_psi;

num    = (P03_psf ./ P0i_psf).^((gamma - 1) / gamma) - 1;
den    = (T03_R ./ T0i_R) - 1;
eta_c  = num ./ den;

figure; grid on; box on; hold on;
plot(mdot, eta_c, 's-', 'LineWidth', 1.5, 'MarkerSize', 7);
xlabel('Air Flow');
ylabel('Compressor Efficiency');
title('Compressor Efficiency vs Air Flow');
text(mdot, eta_c, compose(' %d%%', pctRPM), 'VerticalAlignment','bottom');

% ============================================================
% THRUST CALCULATIONS
% ============================================================
GPS            = [0.00341 0.00383 0.00379 0.00410 0.00433 0.00441 0.00452]; % fuel flow (gal/s)
mdot_f_slugps  = GPS .* rho_f_lbgal .* slug_per_lbm;                        % slug/s

T05_F = [633 625 621 640 644 693 757];
T05_R = T05_F + 459.67;

% P07_psig = [0.27 0.359 0.484 0.426 0.397 0.349 0.344];
% P07_psia = P07_psig + Pamb_psia;

w_open_in = 8 + 1/8;
h_open_in = 5 + 3/8;
A_exit_open = (w_open_in / 12) * (h_open_in / 12); % ft^2

% Nozzle solver for converging nozzle (isentropic to Pe = Pa unless choked)
nozzleThrust = @(m_slugps, P0_psia, T0_R, Aexit_ft2) ...
    local_converging_nozzle_thrust(m_slugps, P0_psia, T0_R, Aexit_ft2, ...
                                   Pa_psia, Pa_psf, gamma_g, R_g, P07_psia);

% ============================================================
% THRUST LOOP
% ============================================================
mdot_exhaust = mdot + mdot_f_slugps;
T_lbf        = zeros(size(pctRPM));

for k = 1:numel(pctRPM)
    T_lbf(k) = thrust(gamma_g, R_hot, A_exit_open, P0i_psf(k), Pi_psf(k), P07_psf(k), P7_psf(k), T0i_R(k), T05_R(k), Pamb_psf);
    % nozzleThrust(mdot_exhaust(k), P07_psia(k), T05_R(k), A_exit_open);
end

Tcorr_lbf = T_lbf .* (T_std_R ./ Ta_R) .* (Pa_psia ./ P_std_psia);

% ============================================================
% THRUST PLOT
% ============================================================
figure; grid on; box on; hold on;
plot(pctRPM, Tcorr_lbf, 'o-', 'LineWidth', 1.5, 'MarkerSize', 7);
xlabel('Percent Throttle (% RPM)');
ylabel('Corrected Thrust, T_c (lbf)');
title('Corrected Thrust vs Percent Throttle (Open Nozzle)');

%% Plots 4-5
% d) EGT (R) vs. engine pressure ratio (EPR)
% e) Corrected thrust (lbf) vs. nozzle setting at 100% throttle

data = readtable("All Sensor Data Export - 2025-10-22 - 16_27 - 16_50.csv");
%T05_R_all = data.TempSensorTo5DegreeF + 459.67;

% EPR = (data.PressureSensorP08PitotImpactPSIg + Pamb_psia) ./...
%     (data.PressureSensorP01InletPitotImpactPSIg + Pamb_psia);
EPR = P07_psia ./ P0i_psia;

% ============================================================
% EGT PLOT
% ============================================================
figure; grid on; box on; hold on;
plot(EPR, T05_R, 'o-', 'LineWidth', 1.5, 'MarkerSize', 7);
xlabel('EPR');
ylabel('EGT (R)');
title('EGT (R) vs. Engine Pressure Ratio (EPR)');

nozzleA = [0.098, 0.134, A_exit_open];
for i = 1:numel(nozzleA)
    k = (7);
    Tn_lbf(i) = thrust(gamma_g, R_hot, nozzleA(i), P0i_psf(k), Pi_psf(k), P07_psf(k), P7_psf(k), T0i_R(k), T05_R(k), Pamb_psf);
end


figure; grid on; box on; hold on;
plot(nozzleA, Tn_lbf, 'o-', 'LineWidth', 1.5, 'MarkerSize', 7);
xlabel('Nozzle');
ylabel('Thrust');
title('Corrected Thrust vs. Nozzle Setting (100% Thrust)');

%% Functions

function T_lbf = thrust(gamma, R, Ae, P0i_psf, Pi_psf, P0e_psf, Pe_psf, T0i_R, T0e_R, Pamb_psf)
%T_lbf = thrust(gamma, R, Ae, P0i_psf, Pi_psf, P0e_psf, Pe_psf, T0i_R, T0e_R)
% Uses rankine for temp, and psf for pressure
% e is exit
% i is inlet
    A1 = pi*(5.5/2/12)^2; % ft^2 %hardcoding inlet area, change this to use function


    M1 = sqrt( (2/(gamma - 1)) * ((P0i_psf./Pi_psf).^((gamma - 1)/gamma) - 1) );
    M8 = sqrt( (2/(gamma - 1)) * ((P0e_psf./Pe_psf).^((gamma - 1)/gamma) - 1) );
    
    mdoti = (A1 * P0i_psf ./ sqrt(T0i_R)) * sqrt(gamma / R) .* ...
           M1 .* (1 + (gamma - 1)./2 .* M1.^2).^(- (gamma + 1) / (2 * (gamma - 1)));
    
    mdote = (Ae * P0e_psf ./ sqrt(T0e_R)) * sqrt(gamma / R) .* ...
           M8 .* (1 + (gamma - 1)./2 .* M8.^2).^(- (gamma + 1) / (2 * (gamma - 1)));
    
    Vi = M1.*sqrt(gamma*R*T0i_R);
    Ve = M8.*sqrt(gamma*R*T0e_R);
    
    T_lbf = mdote.*Ve - mdoti.*Vi - (Pe_psf - Pamb_psf)*Ae;
end


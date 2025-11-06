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
R_cold          = 1716; % ft*lbf/(slug*R)
R_hot           = 1544; % ft*lbf/(slug*R)

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
P7_psia  = P7_psig  + Pamb_psia;

Pamb_psf = Pamb_psia * psf_per_psi;
P7_psf   = P7_psia   * psf_per_psi;
P07_psf  = P07_psia  * psf_per_psi;
Pi_psf   = Pi_psia   * psf_per_psi;      % psf
P0i_psf  = P0i_psia  * psf_per_psi;      % psf
P03_psf  = P03_psia  * psf_per_psi;      % psf

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
gamma_g = 1.33; % hot exhaust (for completeness; nozzle uses R_hot below)
R_g     = 1716;

% Ambient for thrust corrections
Ta_R        = Ta_F + 459.67;

% Jet-A fuel density
rho_f_lbgal   = 6.7;
g_c           = 32.174;
slug_per_lbm  = 1 / g_c;

% Temperatures
T03_F = [197 247 282 314 342 366 396];
T03_R = T03_F + 459.67;

num    = (P03_psf ./ P0i_psf).^((gamma - 1) / gamma) - 1;
den    = (T03_R ./ T0i_R) - 1;
eta_c  = num ./ den;

figure; grid on; box on; hold on;
plot(mdot, eta_c, 's-', 'LineWidth', 1.5, 'MarkerSize', 7);
xlabel('Air Flow (slug/s)');
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

% Nozzle opening geometry (largest/open setting from your rig)
w_open_in = 8 + 1/8;
h_open_in = 5 + 3/8;
A_exit_open = (w_open_in / 12) * (h_open_in / 12); % ft^2

% ============================================================
% THRUST LOOP (uses consistent nozzle solution)
% ============================================================
T_lbf = zeros(size(pctRPM));
for k = 1:numel(pctRPM)
    m_in   = mdot(k);
    m_fuel = mdot_f_slugps(k);
    m_out  = m_in + m_fuel;

    % inlet velocity should use static temperature
    V_in_k = M_in(k) * sqrt(gamma * R_cold * Ti_R(k));

    % Solve nozzle given area, total exit conditions, and mass flow
    [Ve_k, Pe_psf_k] = solve_nozzle_exit_from_area( ...
        A_exit_open, P07_psf(k), T05_R(k), 1.33, R_hot, m_out, Pamb_psf);

    % Momentum + pressure thrust (correct sign on pressure term)
    T_lbf(k) = m_out * Ve_k - m_in * V_in_k + (Pe_psf_k - Pamb_psf) * A_exit_open;
end

% Standard-day correction
Tcorr_lbf = T_lbf .* (T_std_R ./ Ta_R) .* (Pamb_psia ./ P_std_psia);

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
EPR = (data.PressureSensorP08PitotImpactPSIg + Pamb_psia) ./ (data.PressureSensorP01InletPitotImpactPSIg + Pamb_psia);

% data manipulation ------------------------------------------

% === Upward-curving regression on EGT vs EPR ===
x = data.TempSensorTo5DegreeF + 459.67;   % EGT (R)
y = EPR;                                   % EPR

% keep finite points only
good = isfinite(x) & isfinite(y);
x = x(good); y = y(good);

% --- 1) Bin + median to denoise (robust to outliers)
nbins  = 60;                        % tweak if you want smoother/edgier
edges  = linspace(min(x), 1400, nbins+1);
[~,~,bin] = histcounts(x, edges);

% median EPR in each EGT bin (and mean EGT for x-position)
y_med = accumarray(bin(bin>0), y(bin>0), [], @median, NaN);
x_avg = accumarray(bin(bin>0), x(bin>0), [], @mean,   NaN);

mask  = ~isnan(x_avg) & ~isnan(y_med);
x_fit = x_avg(mask);
y_fit = y_med(mask);

% --- 2) Quadratic fit (tends to give a gentle "curves up" shape)
p2      = polyfit(x_fit, y_fit, 5);     % y ~ p2(1)x^2 + p2(2)x + p2(3)
x_line  = linspace(min(x_fit), max(x_fit), 300);
y_line  = polyval(p2, x_line);

% If curvature accidentally ends up negative (rare after binning), refit cubic
if p2(1) < 0
    p3 = polyfit(x_fit, y_fit, 3);
    y_line = polyval(p3, x_line);
end

% --- 3) Plot the denoised medians and the regression curve


% ------------------------------------------------------------

% ============================================================
% EGT PLOT
% ============================================================
figure; grid on; box on; hold on;
plot(data.TempSensorTo5DegreeF + 459.67, EPR, '.', 'LineWidth', 1.5, 'MarkerSize', 7);
plot(x_fit, y_fit, 'ko', 'MarkerSize', 4, 'HandleVisibility', 'off'); % binned medians
plot(x_line, y_line, 'LineWidth', 2, 'DisplayName', 'Quadratic fit');
ylabel('EPR');
xlabel('EGT (R)');
title('EGT (R) vs. Engine Pressure Ratio (EPR)');
legend('Raw Data', 'Polynomial Fit')

% ============================================================
% Corrected Thrust vs Nozzle Setting (100% throttle)
% ============================================================
nozzleA = [0.098, 0.134, A_exit_open];

k = 7; % 100% RPM row
m_in_7   = mdot(k);
m_fuel_7 = mdot_f_slugps(k);
m_out_7  = m_in_7 + m_fuel_7;

V_in_7 = M_in(k) * sqrt(gamma * R_cold * Ti_R(k));

Tn_lbf = zeros(size(nozzleA));
for i = 1:numel(nozzleA)
    Ae = nozzleA(i);
    [Ve_i, Pe_psf_i] = solve_nozzle_exit_from_area( ...
        Ae, P07_psf(k), T05_R(k), 1.33, R_hot, m_out_7, Pamb_psf);

    Tn_lbf(i) = m_out_7 * Ve_i - m_in_7 * V_in_7 + (Pe_psf_i - Pamb_psf) * Ae;
end

figure; grid on; box on; hold on;
plot(nozzleA, Tn_lbf, 'o-', 'LineWidth', 1.5, 'MarkerSize', 7);
xlabel('Nozzle Exit Area (ft^2)');
ylabel('Thrust (lbf)');
title('Corrected Thrust vs. Nozzle Setting (100% Throttle)');

%% Functions

function [Ve, Pe_psf, M8] = solve_nozzle_exit_from_area(Ae, P0e_psf, T0e_R, gamma, R, mdot_req, Pamb_psf)
% Returns exit velocity Ve (ft/s), exit static pressure Pe_psf (psf), and exit Mach M8
% for a converging nozzle with given area Ae, total P0e, total T0e and required mass flow.
% Enforces choking (M=1) if required mass flow exceeds choked capacity.

    % Mass-flow function G(M)
    G = @(M) M .* (1 + (gamma-1)/2 .* M.^2) .^ (-(gamma+1)/(2*(gamma-1)));

    % Scaling constant for mdot
    K = (Ae * P0e_psf / sqrt(T0e_R)) * sqrt(gamma / R);

    % Choked capacity at M=1
    mdot_star = K * G(1.0);

    if mdot_req >= mdot_star * 0.999 % treat as choked
        M8 = 1.0;
        Te_R  = T0e_R / (1 + (gamma-1)/2 * M8^2);      % = T0e * 2/(γ+1)
        Pe_psf = P0e_psf * (Te_R / T0e_R)^(gamma/(gamma-1));
    else
        % Unchoked: find subsonic M in (0,1) such that K*G(M)=mdot_req
        target = mdot_req / K;
        f = @(M) G(M) - target;
        % robust bracket solve
        M8 = fzero(f, [1e-6, 0.999]);
        Te_R  = T0e_R / (1 + (gamma-1)/2 * M8^2);
        Pe_psf = P0e_psf * (Te_R / T0e_R)^(gamma/(gamma-1));
    end

    % Exit velocity
    a_e = sqrt(gamma * R * Te_R);
    Ve  = M8 * a_e;

    % For a simple static test assumption: clamp to ambient if Pe would underexpand
    if Pe_psf < Pamb_psf
        Pe_psf = Pamb_psf;
    end
end

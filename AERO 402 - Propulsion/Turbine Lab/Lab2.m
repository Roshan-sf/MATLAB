%% Lab 2
clc
close all
clear
%%
%data = readmatrix('All Sensor Data Export - 2025-10-22 - 16_27 - 16_50.csv.xlsx');
gamma = 1.4;
R = 1716; 
A1 = pi*(5.5/2/12)^2; % ft^2
w_open_in = 8 + 1/8;
h_open_in = 5 + 3/8;
A8 = (w_open_in / 12) * (h_open_in / 12); % ft^2;
% 
% Pamb = data(:, 3)*144; % PSFA
% To1F = data(:, 14); % F
% To1 = To1F + 459.67; % R
% P1g = data(:, 11)*144; % PSFG
% P1 = P1g + Pamb; % PSFA
% Po1g = data(:, 8)*144; % PSFG
% Po1 = Po1g + Pamb; % PSFA
% LoadCell = data(:, 7); % lbf
% 
% To5F = data(:, 17); % F
% To5 = To5F + 459.67; % R
% Po8g = data(:, 10)*144; % PSFG
% Po8 = Po8g + Pamb; % PSFA
% P8g = data(:, 12)*144; % PSFG
% P8 = P8g + Pamb; % PSFA

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
T05_F = [633 625 621 640 644 693 757];
T05_R = T05_F + 459.67;

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


M1 = sqrt( (2/(gamma - 1)) * ((P0i_psf./Pi_psf).^((gamma - 1)/gamma) - 1) );
M8 = sqrt( (2/(gamma - 1)) * ((P07_psf./P7_psf).^((gamma - 1)/gamma) - 1) );

mdot1 = (A1 * P0i_psf ./ sqrt(T0i_R)) * sqrt(gamma / R) .* ...
       M1 .* (1 + (gamma - 1)./2 .* M1.^2).^(- (gamma + 1) / (2 * (gamma - 1)));

mdot8 = (A8 * P07_psf ./ sqrt(T05_R)) * sqrt(gamma / R) .* ...
       M8 .* (1 + (gamma - 1)./2 .* M8.^2).^(- (gamma + 1) / (2 * (gamma - 1)));

V1 = M1.*sqrt(gamma*R*T0i_R);
V8 = M8.*sqrt(gamma*R*T05_R);

F = mdot8.*V8 - mdot1.*V1 - (P7_psf - Pamb_psf)*A8;

figure
plot(F)
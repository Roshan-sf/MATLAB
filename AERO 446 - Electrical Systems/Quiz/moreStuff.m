clc; clear;

% ------------------------------
% Given Parameters
% ------------------------------
fc = 4e9;                  % Carrier frequency (Hz)
BW = 50e6;                 % Bandwidth (Hz)
P_TX_Watts = 75;          % Transmit power in watts
OBO_dB = 3;               % Output Back-Off (dB)
Tx_line_loss_dB = 2;      % Transmit line loss (dB)

G_RX_dB = 21;             % Receive antenna gain (dB)
NF_dB = 3;                % Noise figure of receiver (dB)
C_N_required_dB = 10;     % Required C/N (dB)

R_geo = 42164.18e3;       % GEO orbital radius (m)
lambda = 3e8 / fc;        % Wavelength (m)

% ------------------------------
% Transmit Power after OBO and losses
% ------------------------------
P_out_dBW = 10*log10(P_TX_Watts) - OBO_dB;  % Convert to dBW and subtract OBO
P_out_dBW = P_out_dBW - Tx_line_loss_dB;   % Apply transmitter line loss

% ------------------------------
% Atmospheric and Path Losses (from lecture slide)
% ------------------------------
L_fs_dB = 196.39;          % Free space loss
L_atm_dB = 7.45;           % Atmospheric absorption
L_pol_dB = 0.3;            % Polarization loss
L_radome_dB = 1;           % Radome loss
L_total_dB = L_fs_dB + L_atm_dB + L_pol_dB + L_radome_dB;

% ------------------------------
% Receiver G/T
% ------------------------------
T_s = 290 * 10^(NF_dB/10);         % System noise temperature (K)
G_T_dB_per_K = G_RX_dB - 10*log10(T_s);  % G/T in dB/K

% ------------------------------
% Compute Required EIRP
% ------------------------------
C_N_dB = C_N_required_dB;
EIRP_required_dB = C_N_dB - G_T_dB_per_K + L_total_dB - 228.6 + 10*log10(BW);

% ------------------------------
% Required Transmit Antenna Gain
% ------------------------------
G_T_dBi = EIRP_required_dB - P_out_dBW;

% ------------------------------
% Solve for Antenna Diameter
% G = 10*log10(η * (pi*D/lambda)^2)
% => D = (lambda/pi) * sqrt(10^(G/10) / η)
% ------------------------------
eta = 0.55;  % Efficiency
D = (lambda/pi) * sqrt(10^(G_T_dBi/10) / eta);

% ------------------------------
% Display Results
% ------------------------------
fprintf('--- Earth Receive Antenna Sizing ---\n');
fprintf('Required Transmit EIRP: %.2f dBW\n', EIRP_required_dB);
fprintf('Required Transmit Antenna Gain: %.2f dBi\n', G_T_dBi);
fprintf('Required Earth Receive Antenna Diameter: %.2f meters\n', D);

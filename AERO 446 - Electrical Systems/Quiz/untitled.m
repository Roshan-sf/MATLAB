% Constants
c = 3e8;                       % Speed of light (m/s)
k = 1.38e-23;                  % Boltzmann constant (J/K)
BW = 50e6;                     % Bandwidth in Hz
f = 4e9;                       % Frequency in Hz
lambda = c / f;               % Wavelength in meters

% Transmit Power and Losses
P_tx_W = 75;                  % Transmit Power in Watts
OBO_dB = 3;                   % Output Back-Off in dB
P_tx_dBW = 10*log10(P_tx_W) - OBO_dB;  % Effective power in dBW
L_tx = 2;                     % Transmitter loss in dB

% EIRP
G_tx = 0;                     % Placeholder (antenna gain unknown yet)
EIRP = P_tx_dBW - L_tx + G_tx;  % Will update later

% Receiver G/T Calculation
G_r_dB = 21;                  % Receiver gain in dB
F_dB = 3;                     % Noise figure in dB
T_sys = 290 * 10^(F_dB/10);   % System noise temperature in Kelvin
G_over_T = G_r_dB - 10*log10(T_sys);  % G/T in dB

% Free space path loss (FSPL)
d_km = sqrt((42164.18)^2 + (6378.18)^2 - 2*42164.18*6378.18*cosd(60));
d_m = d_km * 1e3;
FSPL_dB = 20*log10(4*pi*d_m / lambda);  % Free-space path loss in dB

% Link Budget: C/N
C_N = EIRP + G_over_T - FSPL_dB + 228.6 - 10*log10(BW);

% Antenna Gain Requirement
C_N_required = 10;  % dB
margin = 3;         % dB margin
C_N_target = C_N_required + margin;

% Backsolve for required G_tx
G_tx_required = C_N_target - (P_tx_dBW - L_tx + G_over_T - FSPL_dB + 228.6 - 10*log10(BW));
fprintf('Required Transmit Antenna Gain: %.2f dB\n', G_tx_required);

% Parabolic antenna diameter
efficiency = 0.55;  % Typical efficiency
G_tx_linear = 10^(G_tx_required/10);
D = sqrt((G_tx_linear * lambda^2) / (pi^2 * efficiency));
fprintf('Required Antenna Diameter: %.2f meters\n', D);

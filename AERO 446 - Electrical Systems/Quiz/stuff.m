%% Battery Calculations

% -------------------------------
% Variables (example placeholder values)

P_out = 500;            % Power needed from battery [W]
t_discharge = 2;        % Duration battery is used [hr]
DOD = 0.8;              % Depth of Discharge [fraction]
eta = 0.9;              % Load Transmission Efficiency [fraction]

t_charge = 3;           % Duration battery is charged [hr]
V_charge = 4.2;         % Charge voltage [V]
V_discharge = 3.6;      % Discharge voltage [V]

V_tn = 3.8;             % Thermoneutral voltage [V]

% -------------------------------
% 1. Battery Capacity [Wh]
E = (P_out * t_discharge) / (DOD * eta);  

% -------------------------------
% 2. Power needed to charge battery [W]
P_in = P_out * (t_discharge * V_charge) / (t_charge * V_discharge);

% -------------------------------
% 3. Battery Heat Dissipation [W]
P_thermal = (P_out * V_tn / V_discharge) - P_out;

% -------------------------------
% 4. Battery Efficiency [fraction]
eta_eff = V_charge / V_discharge;

% Display results
fprintf('Battery Capacity: %.2f Wh\n', E);
fprintf('Power needed to charge battery: %.2f W\n', P_in);
fprintf('Heat dissipation: %.2f W\n', P_thermal);
fprintf('Battery Efficiency: %.2f %%\n', eta_eff * 100);

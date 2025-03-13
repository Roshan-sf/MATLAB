% Steady-State Heat Transfer Calculations for Metal Rods
clc; 
clear; 
close all;

% Given Data
D_rod = 1/4 / 12; % Rod diameter in feet (converted from inches)
T_water = 200; % Water temperature (F)
T_air = 130; % Ambient air temperature (F)
g = 32.174; % Gravity (ft/s^2)
Pr = 0.71; % Prandtl number for air
nu = 1.6e-4; % Kinematic viscosity of air (ft^2/s)

% Material Properties (Thermal Conductivity in Btu/hr·ft·°F)
materials = {'Aluminum', 'Copper', 'Carbon Steel'};
k_values = [118, 223, 26]; % Thermal conductivity (Btu/hr·ft·°F)

% Mass in grams (g) for each material
masses_g = [26.6, 88.3, 77.9]; % Given mass values in grams
masses = masses_g * 0.00220462; % Convert mass to lbm

% Rod Length for Each Material
L_rods = [9.5/12, 9.5/12, 9.5/12]; % Rod length for each material (ft)

% Loop through each material
for i = 1:length(materials)
    L_rod = L_rods(i); % Get rod length for this material
    k = k_values(i); % Get thermal conductivity for this material
    mass = masses(i); % Get converted mass for this material
    
    % Compute common parameters
    A_cross = pi * (D_rod^2) / 4; % Cross-sectional area for conduction (ft^2)
    A_surface = pi * D_rod * L_rod; % Surface area for convection (ft^2)
    L_c = L_rod; % Characteristic length for natural convection
    beta = 1 / (T_air + 459.67); % Thermal expansion coefficient (1/R) using Rankine conversion

    % Compute Grashof Number (Gr)
    Gr = (g * beta * (T_water - T_air) * (L_c^3)) / (nu^2);

    % Compute Nusselt Number (Nu)
    Nu = 0.59 * (Gr * Pr)^(1/4);

    % Compute Convective Heat Transfer Coefficient (h)
    h = (Nu * k) / L_c;

    % Compute Thermal Resistances
    R_conduction = L_rod / (k * A_cross); % Conduction resistance
    R_convection = 1 / (2 * pi * (D_rod / 2) * h); % Convection resistance

    % Compute Heat Transfer Rate (q)
    q = (T_water - T_air) / (R_conduction + R_convection);

 % Display Results
    fprintf('Material: %s\n', materials{i});
    fprintf('---------------------------------------------\n');
    fprintf('| Rod Length (L_rod)         : %.3f ft\n', L_rods(i));
    fprintf('| Mass (m)                   : %.3f lbm\n', masses(i));
    fprintf('| Convective Heat Transfer Coeff. (h) : %.3f Btu/hr·ft²·°F\n', h);
    fprintf('| Thermal Resistance (R_cond) : %.5f hr·ft²·°F/Btu\n', R_conduction);
    fprintf('| Thermal Resistance (R_conv) : %.5f hr·ft²·°F/Btu\n', R_convection);
    fprintf('| Total Thermal Resistance (R_total) : %.5f hr·ft²·°F/Btu\n', R_conduction + R_convection);
    fprintf('| Heat Transfer Rate (q)      : %.3f Btu/hr\n', q);
    fprintf('---------------------------------------------\n\n');
end

clear

% Given Data
T_c = 3450; % Combustion chamber temperature (Rankine)
P_c = 110; % Combustion chamber pressure (psi)
T_inf = 130 + 459.67; % Ambient temperature (Rankine)
t_burn = 1.6; % Burn time (seconds)
t_cool = 82; % Cooling time (seconds)
t_total = t_burn + t_cool; % Total time

% Thermal properties
gamma = 1.25; % Specific heat ratio
k_g = 0.06332; % Thermal conductivity of combustion gas (Btu/hr·ft·°F)
k_p = 0.15; % Thermal conductivity of phenolic casing (Btu/hr·ft·°F)
D_outer = 0.69 / 12; % Outer diameter of motor casing (feet)
t_p = 1/16 / 12; % Casing thickness (feet)
D_inner = D_outer - 2*t_p; % Inner diameter of casing (feet)
L_c = 7 / 12; % Characteristic length (feet)
Pr = 0.7496; % Prandtl number
Re = 50000; % Approximate Reynolds number

% Step 2: Compute Nozzle Throat Conditions
T_star = T_c / (1 + (gamma - 1) / 2); % Nozzle throat temperature
P_star = P_c * (1 / (1 + (gamma - 1) / 2))^(gamma / (gamma - 1)); % Nozzle throat pressure

% Step 3: Compute Convective Heat Transfer Coefficient (h)
Nu = 0.023 * (Re^0.8) * (Pr^0.4); % Nusselt number
h = (Nu * k_g) / D_inner; % Convective heat transfer coefficient

% Step 4: Compute Thermal Resistance
A = pi * D_outer * L_c; % Surface area for heat transfer
R_cond = t_p / (k_p * A); % Conduction resistance
R_conv = 1 / (h * A); % Convection resistance
R_total = R_cond + R_conv; % Total thermal resistance

% Step 5: Compute Heat Transfer Rate During Burn (q_burn)
q_burn = (T_c - T_inf) / R_total; % Heat transfer rate during burn

% Step 6: Compute Heat Transfer Rate During Cooling (q_cool)
q_cool = q_burn / 5; % Cooling heat transfer rate (assumed to be 1/5 of q_burn)

% Step 7: Compute Time-Averaged Heat Transfer Rate (q_avg)
q_avg = (q_burn * (t_burn / t_total)) + (q_cool * (1 - (t_burn / t_total)));

% Apply Correction Factor (0.9) to q_avg
q_avg = q_avg * 0.9;

% Step 8: Compute Outer Surface Temperature (T3)
T3 = T_c - q_avg * R_total;

% Display Results in a Readable Format
fprintf('    Heat Transfer Analysis for Rocket Motor  \n');

fprintf('Nozzle Throat Temperature (T*): %.2f Rankine\n', T_star);
fprintf('Nozzle Throat Pressure (P*): %.2f psi\n', P_star);

fprintf('Nusselt Number (Nu): %.3f\n', Nu);
fprintf('Convective Heat Transfer Coefficient (h): %.3f Btu/hr·ft²·°F\n', h);

fprintf('Thermal Resistance Values:\n');
fprintf('  - Conduction Resistance (R_cond): %.5f hr·ft²·°F/Btu\n', R_cond);
fprintf('  - Convection Resistance (R_conv): %.5f hr·ft²·°F/Btu\n', R_conv);
fprintf('  - Total Thermal Resistance (R_total): %.5f hr·ft²·°F/Btu\n', R_total);

fprintf('Heat Transfer Rates:\n');
fprintf('  - During Burn (q_burn): %.3f Btu/hr\n', q_burn);
fprintf('  - During Cooling (q_cool): %.3f Btu/hr\n', q_cool);
fprintf('  - Time-Averaged (q_avg): %.3f Btu/hr\n', q_avg);

fprintf('Final Outer Surface Temperature (T3): %.2f Rankine\n', T3);

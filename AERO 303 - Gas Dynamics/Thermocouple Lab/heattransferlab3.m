% Steady-State Heat Transfer Calculations for Metal Rods (Metric Units)
clc; clear; close all;

% Given Data (Metric)
D_rod = (1/4) * 0.0254; % Rod diameter in meters (converted from inches)
T_water = (200 - 32) * 5/9 + 273.15; % Water temperature in Kelvin
T_air = (130 - 32) * 5/9 + 273.15; % Ambient air temperature in Kelvin
g = 9.81; % Gravity in m/s²
Pr = 0.71; % Prandtl number for air
nu = 18.9e-6; % Kinematic viscosity in m²/s

% Material Properties (Thermal Conductivity in W/m·K)
materials = {'Aluminum', 'Copper', 'Carbon Steel'};
k_values = [118, 223, 26] * 1.731; % Convert Btu/hr·ft·°F to W/m·K

% Mass in grams (g) converted to kilograms (kg)
masses_g = [26.6, 88.3, 77.9]; % Given mass values in grams
masses = masses_g * 0.001; % Convert mass to kg

% Rod Length for Each Material (feet to meters)
L_rods = [9.5/12, 9.5/12, 9.5/12] * 0.3048; % Convert to meters

% Display header once
fprintf('=============================================\n');
fprintf('       Heat Transfer Analysis Results (SI)  \n');
fprintf('=============================================\n\n');

% Loop through each material
for i = 1:length(materials)
    L_rod = L_rods(i); % Get rod length for this material (m)
    k = k_values(i); % Get thermal conductivity for this material (W/m·K)
    mass = masses(i); % Get converted mass for this material (kg)
    
    % Compute common parameters
    A_cross = pi * (D_rod^2) / 4; % Cross-sectional area for conduction (m²)
    A_surface = pi * D_rod * L_rod; % Surface area for convection (m²)
    L_c = L_rod; % Characteristic length for natural convection (m)
    beta = 1 / T_air; % Thermal expansion coefficient (1/K) (Use Kelvin)

    % Compute Grashof Number (Gr) using Kelvin temperatures
    Gr = (g * beta * (T_water - T_air) * (L_c^3)) / (nu^2);

    % Compute Nusselt Number (Nu)
    Nu = 0.59 * (Gr * Pr)^(1/4);

    % Compute Convective Heat Transfer Coefficient (h)
    h = (Nu * 0.0262) / L_c; % W/m²·K

    % Compute Thermal Resistances
    R_conduction = L_rod / (k * A_cross); % Conduction resistance (m²·K/W)
    R_convection = 1 / (h * A_surface); % Convection resistance (m²·K/W)
    R_total = R_conduction + R_convection; % Total resistance (m²·K/W)

    % Compute Heat Transfer Rate (q)
    q = (T_water - T_air) / R_total; % Heat transfer rate (W)

    % Display Results for the current material
    fprintf('Material: %s\n', materials{i});
    fprintf('---------------------------------------------\n');
    fprintf('| Rod Length (L_rod)         : %.3f m\n', L_rod);
    fprintf('| Mass (m)                   : %.3f kg\n', mass);
    fprintf('| Convective Heat Transfer Coeff. (h) : %.3f W/m²·K\n', h);
    fprintf('| Thermal Resistance (R_cond) : %.5f m²·K/W\n', R_conduction);
    fprintf('| Thermal Resistance (R_conv) : %.5f m²·K/W\n', R_convection);
    fprintf('| Total Thermal Resistance (R_total) : %.5f m²·K/W\n', R_total);
    fprintf('| Heat Transfer Rate (q)      : %.3f W\n', q);
    fprintf('---------------------------------------------\n\n');
end

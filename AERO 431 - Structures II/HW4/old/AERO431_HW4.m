% AERO 431
% HW #4 
% Sebastian Macias-Gonzalez
% Due June 1st, 2025

clear; clc; close all;

%% Question 1
disp("-----")
disp("Question 1")

% Constants
E = 70e9;             % Young's Modulus [Pa]
G = 26e9;             % Shear Modulus [Pa]
sigma_y = 280e6;      % Yield strength [Pa]
t = 0.002;            % Skin thickness [m]
T = 0.005;            % Flange thickness [m]
h = 0.1;              % Web height [m]
w = 0.05;             % Flange width [m]
nu = E / (2 * G) - 1; % Poisson's ratio

% Define L and b ranges
L_vals = linspace(0.1, 0.5, 6);
b_vals = linspace(0.1, 0.3, 6);

% Storage
results = [];

% Buckling coefficient function (for simply supported plates)
k_plate = @(aspect, m) (pi^2) * (m^2 + (aspect * m)^2);

% Plate buckling stress function
buckling_stress = @(E, k, t_val, b_val) (k * pi^2 * E / (12 * (1 - nu^2))) * (t_val / b_val)^2;

% Loop through all combinations
for i = 1:length(L_vals)
    for j = 1:length(b_vals)
        L = L_vals(i);
        b = b_vals(j);
        m = 1;

        % Aspect ratios
        aspect_skin = L / b;
        aspect_web = h / b;
        aspect_flange = L / (w / 2);
        aspect_Iplate = h / b;

        % Buckling coefficients
        k_skin = k_plate(aspect_skin, m);
        k_web = k_plate(aspect_web, m);
        k_flange = k_plate(aspect_flange, m);
        k_Iplate = k_plate(aspect_Iplate, m);

        % Buckling stresses (plate theory)
        sigma_skin = min(buckling_stress(E, k_skin, t, b), sigma_y);
        sigma_web = min(buckling_stress(E, k_web, t, b), sigma_y);
        sigma_flange = min(buckling_stress(E, k_flange, T, w/2), sigma_y);
        sigma_Iplate = min(buckling_stress(E, k_Iplate, T, b), sigma_y);

        % I-beam buckling (Euler)
        I_beam = 2 * (w * T^3) / 12 + 2 * w * T * (h/2)^2;
        A_beam = 2 * w * T + h * t;
        sigma_Euler = min((pi^2 * E * I_beam) / (A_beam * L^2), sigma_y);

        % Critical stress for each method
        sigmas_plate = [sigma_skin, sigma_web, sigma_flange, sigma_Iplate];
        sigmas_Ibeam = [sigma_skin, sigma_web, sigma_flange, sigma_Euler];

        [sigma_crit_plate, idx_plate] = min(sigmas_plate);
        [sigma_crit_Ibeam, idx_Ibeam] = min(sigmas_Ibeam);

        modes = {'Skin Panel', 'Spar Web', 'Spar Flange', 'I-beam Plate'};

        if idx_Ibeam == 4
            modes{4} = 'I-beam Beam';
        end

        % Store results
        results = [results;
            L, b, ...
            sigma_skin, sigma_web, sigma_flange, sigma_Iplate, sigma_Euler, ...
            sigma_crit_plate, sigma_crit_Ibeam, ...
            idx_plate, idx_Ibeam];
    end
end

% Display as table
headers = {'L (m)', 'b (m)', 'σ_skin (Pa)', 'σ_web (Pa)', 'σ_flange (Pa)', ...
           'σ_Iplate (Pa)', 'σ_Euler (Pa)', 'σ_crit_plate (Pa)', ...
           'σ_crit_Ibeam (Pa)', 'Mode_Plate_Idx', 'Mode_Ibeam_Idx'};

result_table = array2table(results, 'VariableNames', headers);

% Add readable mode names
mode_names = {'Skin Panel', 'Spar Web', 'Spar Flange', 'I-beam Plate/Beam'};
result_table.Mode_Plate = mode_names(result_table.Mode_Plate_Idx)';
result_table.Mode_Ibeam = mode_names(result_table.Mode_Ibeam_Idx)';

% Show final table
disp(result_table)

% Explantion 
% The results show that skin buckling is the dominant failure mode, highlighting
% the importance of proper panel sizing in aircraft structural design. While a
% simple plate model is adequate for analysing short spars, an I-beam model becomes
% necessary as the spar span increases. As the spacing between spars increases (with
% increasing b), the buckling stress decreases—particularly for skin and web panels.
% Similarly, as rib spacing (L) increases, the buckling capacity of flanges and I-beams
% columns also decreases.

%% Question 2
disp("-----")
disp("Question 2")

% Constants
W = 0.5;                % Plate width [m]
sigma = 50e6;           % Applied stress [Pa]
K_IC = 24e6;            % Fracture toughness [Pa*sqrt(m)]
E = 70e9;               % Young's modulus [Pa]

% Paris' law constants
C = 3.15e-11;           % [m/cycle]/(MPa*sqrt(m))^m
m = 3;                  % Paris exponent

% Define beta and K_I functions
beta = @(a) (1.122 - 1.122*(a/W) - 0.820*(a/W).^2 + 3.768*(a/W).^3 - 3.040*(a/W).^4) ./ sqrt(1 - 2*a/W);
K_I = @(a) beta(a) .* sigma .* sqrt(pi .* a);

% Solve for a_max using fzero
f = @(a) K_I(a) - K_IC;
a_max = fzero(f, [1e-6, 0.24]);
disp("Maximum safe crack length a_max = " + num2str(a_max, '%.5f') + " m")

% Euler integration with ΔN = 100 cycles
a0 = a_max / 2;
a = a0;
N = 0;
deltaN = 100;

while a < a_max
    K = K_I(a);
    da_dN = C * (K / 1e6)^m;  % Convert K to MPa*sqrt(m)
    a = a + da_dN * deltaN;
    N = N + deltaN;

    % Avoid overshooting beyond a_max
    if a >= a_max
        break;
    end
end

disp("Cycles to failure with ΔN = 100: " + N + " cycles")

% Euler integration with ΔN = 10 cycles
a = a0;
N = 0;
deltaN = 10;

while a < a_max
    K = K_I(a);
    da_dN = C * (K / 1e6)^m;
    a = a + da_dN * deltaN;
    N = N + deltaN;

    if a >= a_max
        break;
    end
end

disp("Cycles to failure with ΔN = 10: " + N + " cycles")

% Explanation
% Comparing step sizes of 10 and 100 reveals a difference of just 80 cycles
% to failure—only 0.072%. This indicates that increasing the step size has
% a minimal impact on accuracy. Therefore, a step size of 100 provides an
% efficient and sufficiently accurate solution.


% AERO 431 Assignment 4
% 6/1/2025
% Alessandro Tedeschi, Stefan Rosu, Roshan Jaiswal Ferri

%% Workspace prep
clear; 
clc; 
close all;

%% Question 1

% Constants
E = 70e9;             % Young modulus [Pa]
G = 26e9;             % Shear modulus [Pa]
sig_y = 280e6;        % Yield strength [Pa]
nu = (E / (2 * G)) - 1; % Poisson ratio

t = 2e-3;            % Thickness [m]
T = 5e-3;            % Flange thickness [m]
h = 0.1;              % Web height [m]
w = 0.05;             % Flange width [m]

% Initiailize L and b 
L_ = linspace(0.1, 0.5, 6);
b_ = linspace(0.1, 0.3, 6);

% Simply supported plate buckling coefficient 
k_plate = @(AR, m) (pi^2) * (m^2 + (AR * m)^2);

% Simply supported plate buckling stress
buckling_stress = @(E, k, t_, b_) (k * pi^2 * E / (12 * (1 - nu^2))) * (t_ / b_)^2;

% Prepare result array
results = [];

% For loop to iterate through all L and b possibilities
for i = 1:length(L_)

for j = 1:length(b_)

    L = L_(i);
    b = b_(j);
    mode = 1;

    % Calculate aspect ratios for each structural element
    AR_skin = L / b;
    AR_web = h / b;
    AR_flange = L / (w / 2);
    AR_Iplate = h / b;

    % Calculate respective buckling coefficients from aspect ratio
    k_skin = k_plate(AR_skin, mode);
    k_web = k_plate(AR_web, mode);
    k_flange = k_plate(AR_flange, mode);
    k_Iplate = k_plate(AR_Iplate, mode);

    % Calculate buckling stress for each element
    sigma_skin = min(buckling_stress(E, k_skin, t, b), sig_y);
    sigma_web = min(buckling_stress(E, k_web, t, b), sig_y);
    sigma_flange = min(buckling_stress(E, k_flange, T, w/2), sig_y);
    sigma_Iplate = min(buckling_stress(E, k_Iplate, T, b), sig_y);
        I_beam = 2 * (w * T^3) / 12 + 2 * w * T * (h/2)^2;
        A_beam = (2 * w * T) + (h * t);
    sigma_Euler = min((pi^2 * E * I_beam) / (A_beam * L^2), sig_y);

    % Minimum stresses as critical stress value 
    sigmas_plate = [sigma_skin, sigma_web, sigma_flange, sigma_Iplate];
    sigmas_Ibeam = [sigma_skin, sigma_web, sigma_flange, sigma_Euler];

    [sigma_crit_plate, i_plate] = min(sigmas_plate);
    [sigma_crit_Ibeam, i_Ibeam] = min(sigmas_Ibeam);

    components = {'Skin panel', 'Spar web', 'Spar flange', 'I-beam plate'};

    if i_Ibeam == 4
        components{4} = 'I-beam';
    end

    % Store results
    results = [results;
        L, b, ...
        sigma_skin, sigma_web, sigma_flange, sigma_Iplate, sigma_Euler, ...
        sigma_crit_plate, sigma_crit_Ibeam, ...
        i_plate, i_Ibeam];
end

end

% Display as table
headers = {'L [m]', 'b [m]', 'sig_skin [Pa]', 'sig_web [Pa]', 'sig_flange [Pa]', ...
           'sig_Iplate [Pa]', 'sig_Euler [Pa]', 'sig_crit_plate [Pa]', ...
           'sig_crit_Ibeam [Pa]', 'Mode_Plate_Idx', 'Mode_Ibeam_Idx'};

% Organize results
solution = array2table(results, 'VariableNames', headers);
mode_names = {'Skin Panel', 'Spar Web', 'Spar Flange', 'I-beam Plate/Beam'};
solution.Mode_Plate = mode_names(solution.Mode_Plate_Idx)';
solution.Mode_Ibeam = mode_names(solution.Mode_Ibeam_Idx)';

% Display results
disp(solution)










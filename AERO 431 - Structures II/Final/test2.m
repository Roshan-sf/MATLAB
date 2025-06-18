%% Problem 1 - Parts A through E

% Geometry and Load
L = 1.0;              % Length (m)
b = 0.1;              % Width (m)
t_total = 0.025;      % Total thickness (m)
P = 1000;             % Tip load (N)

% Ply Data
Nplies = 5;
t_ply = t_total / Nplies;

% Material Properties (Carbon/Epoxy)
E1 = 140e9;           % Pa
E2 = 10e9;            % Pa
G12 = 7e9;            % Pa
nu12 = 0.3;

%% Part A - Tip Deflection and Twist (Symmetric Layups)
layup_angles = {
    [75 -75 75 -75 75],
    [60 -60 60 -60 60],
    [45 -45 45 -45 45],
    [30 -30 30 -30 30],
    [10 -10 10 -10 10]
};
labels = {'[±75°]', '[±60°]', '[±45°]', '[±30°]', '[±10°]'};

w_tips = zeros(1,5);
twist_tips = zeros(1,5);

disp('----Part A----');
for i = 1:5
    theta = layup_angles{i};
    D = compute_D_matrix(theta, t_ply, E1, E2, G12, nu12);
    D11 = D(1,1); D16 = D(1,3); D66 = D(3,3);
    EI = D11 * b;
    GJ = 4 * D66 * b;
    K = 2 * D16 * b;
    S = (1/12) * D11 * b^3;
    [w_tip, twist_tip] = compute_tip_response(EI, GJ, S, K, L, P);
    w_tips(i) = w_tip;
    twist_tips(i) = twist_tip;
    twist_deg = twist_tip * (180/pi);
    disp([labels{i} ': w_tip = ' num2str(w_tip, '%.4e') ' m, twist = ' num2str(twist_deg, '%.4e') ' deg']);
end

%% Part B - Zero Twist Layup Design
layup = [0 0 0 0 0]; % Found to give zero twist
D = compute_D_matrix(layup, t_ply, E1, E2, G12, nu12);
D11 = D(1,1); D16 = D(1,3); D66 = D(3,3);
EI = D11 * b;
GJ = 4 * D66 * b;
K = 2 * D16 * b;
S = (1/12) * D11 * b^3;
[w_tip, twist_tip] = compute_tip_response(EI, GJ, S, K, L, P);
twist_deg = twist_tip * (180/pi);

disp('----Part B----');
disp(['Custom Layup: ' num2str(layup)]);
disp(['w_tip = ' num2str(w_tip, '%.4e') ' m']);
disp(['twist = ' num2str(twist_deg, '%.4e') ' deg']);
if abs(twist_deg) < 1e-3
    disp('Twist is effectively zero — this layup is a strong candidate.');
else
    disp('Twist is not zero — consider modifying layup or testing others.');
end

%% Part C - Minimize Tip Deflection
layups = {
    [0 15 0 15 0],
    [10 -10 0 -10 10],
    [0 30 0 30 0],
    [0 -45 0 -45 0],
    [0 -60 0 -60 0]
};
labels = {
    '[0 15 0 15 0]',
    '[10 -10 0 -10 10]',
    '[0 30 0 30 0]',
    '[0 -45 0 -45 0]',
    '[0 -60 0 -60 0]'
};

w_tips = zeros(1,5);
twist_tips = zeros(1,5);
disp('----Part C----');

for i = 1:length(layups)
    layup = layups{i};
    label = labels{i};
    D = compute_D_matrix(layup, t_ply, E1, E2, G12, nu12);
    D11 = D(1,1); D16 = D(1,3); D66 = D(3,3);
    EI = D11 * b;
    GJ = 4 * D66 * b;
    K = 2 * D16 * b;
    S = (1/12) * D11 * b^3;
    [w_tip, twist_tip] = compute_tip_response(EI, GJ, S, K, L, P);
    twist_deg = twist_tip * (180/pi);
    w_tips(i) = w_tip;
    twist_tips(i) = twist_deg;

    disp(['Layup ' label ':']);
    disp(['w_tip = ' num2str(w_tip, '%.4e') ' m']);
    disp(['twist = ' num2str(twist_deg, '%.4e') ' deg']);
    disp('');
end

%% Part D - Ply Failure Analysis
S1 = 1500e6; S2 = 40e6; S12 = 70e6;
layup = [0 -60 0 -60 0];
label = '[0 -60 0 -60 0]';

D = compute_D_matrix(layup, t_ply, E1, E2, G12, nu12);
D11 = D(1,1); D16 = D(1,3); D66 = D(3,3);
EI = D11 * b;
GJ = 4 * D66 * b;
K = 2 * D16 * b;
S = (1/12) * D11 * b^3;

load_vals = linspace(50, 5000, 500);
z_mids = linspace(-t_total/2 + t_ply/2, t_total/2 - t_ply/2, Nplies);

failure_load = NaN; failure_ply = NaN; failure_index = NaN;

for P = load_vals
    [a3, b3] = compute_tip_response_raw(EI, GJ, S, K, L, P);
    kappa = [a3; b3; 2*a3];
    
    for i = 1:Nplies
        theta = layup(i);
        z = z_mids(i);
        eps = z * kappa;
        Q = get_Q_matrix(E1, E2, G12, nu12);
        Qbar = transform_Q(Q, theta);
        sigma = Qbar * eps;
        sigma1 = sigma(1); sigma2 = sigma(2); tau12 = sigma(3);

        FI = (sigma1/S1)^2 - (sigma1*sigma2)/(S1^2) + ...
             (sigma2/S2)^2 + (tau12/S12)^2;

        if FI >= 1
            failure_load = P;
            failure_ply = i;
            failure_index = FI;
            break;
        end
    end
    if ~isnan(failure_load)
        break;
    end
end

disp('----Part D----');
disp(['Layup with max twist: ' label]);
if ~isnan(failure_load)
    disp(['Failure occurs at P = ' num2str(failure_load, '%.2f') ' N']);
    disp(['Critical ply index: ' num2str(failure_ply)]);
    disp(['Failure index: ' num2str(failure_index, '%.3f')]);
else
    disp('No failure detected up to maximum load.');
end

%% Part E - Natural Frequency Comparison
rho = 1615; % kg/m^3

layups = {
    [75 -75 75 -75 75],
    [60 -60 60 -60 60],
    [45 -45 45 -45 45],
    [30 -30 30 -30 30],
    [10 -10 10 -10 10],
    [0 0 0 0 0],
    [0 15 0 15 0],
    [10 -10 0 -10 10],
    [0 30 0 30 0],
    [0 -45 0 -45 0],
    [0 -60 0 -60 0]
};
labels = {'[±75°]', '[±60°]', '[±45°]', '[±30°]', '[±10°]', ...
          '[0°]', '[0/15]', '[±10/0]', '[0/30]', '[0/-45]', '[0/-60]'};

frequencies_Hz = zeros(1, length(layups));

disp('----Part E----');
for i = 1:length(layups)
    theta = layups{i};
    D = compute_D_matrix(theta, t_ply, E1, E2, G12, nu12);
    D11 = D(1,1); D66 = D(3,3);
    EI = D11 * b;
    m = rho * b * t_total * L;

    num = EI * integral_phi_dd_squared(L);
    den = m * integral_phi_squared(L);
    omega_rad = sqrt(num / den);
    freq_Hz = omega_rad / (2*pi);
    frequencies_Hz(i) = freq_Hz;

    disp([labels{i} ': f = ' num2str(freq_Hz, '%.2f') ' Hz']);
end

%% Problem 2 - Buckling of Composite Skin Panel

% Material properties (Graphite/Epoxy Composite)
E1 = 181e9;         % Pa
E2 = 10.3e9;        % Pa
G12 = 7.17e9;       % Pa
nu12 = 0.28;

% Geometry
t_total = 0.005;    % Total laminate thickness (m)
Nplies = 5;
t_ply = t_total / Nplies;

% Grid for rib and spar spacing
L_vals = linspace(0.1, 0.5, 6);   % Rib spacing (m)
b_vals = linspace(0.1, 0.3, 6);   % Spar spacing (m)
modes = 1:6;                      % Mode numbers

% Layups to analyze
layups = {
    [0 0 0 0 0],            % [0]_s
    [0 90 90 90 0],         % [0/90]_s
    [90 0 0 0 90]           % [90/0]_s
};
labels = {'[0]_s', '[0/90]_s', '[90/0]_s'};

% Preallocate results
Ncr_all = cell(1,3);

for k = 1:3
    layup = layups{k};
    D = compute_D_matrix(layup, t_ply, E1, E2, G12, nu12);
    D11 = D(1,1); D12 = D(1,2); D22 = D(2,2); D66 = D(3,3);

    Ncr_grid = zeros(length(L_vals), length(b_vals));
    min_stress_vals = zeros(36, 1);
    L_vals_flat = zeros(36, 1);
    b_vals_flat = zeros(36, 1);
    idx_flat = 1;

    for i = 1:length(L_vals)
        for j = 1:length(b_vals)
            L = L_vals(i);
            b = b_vals(j);
            R = L / b;

            % Compute Ncr for each mode
            Ncr_vals = zeros(size(modes));
            for m = modes
                Ncr_vals(m) = (pi^2 / (m^2 * L^2)) * ...
                    (D11 * m^4 + 2 * (D12 + 2 * D66) * m^2 * R^2 + D22 * R^4);
            end

            Ncr_grid(i, j) = min(Ncr_vals);
            min_stress_vals(idx_flat) = Ncr_grid(i,j) / t_total / 1e6; % MPa
            L_vals_flat(idx_flat) = L;
            b_vals_flat(idx_flat) = b;
            idx_flat = idx_flat + 1;
        end
    end

    Ncr_all{k} = Ncr_grid;

    % Contour Plot
    figure;
    contourf(L_vals, b_vals, Ncr_grid' / t_total / 1e6, 20, 'LineColor', 'none');
    colorbar;
    xlabel('Rib spacing L (m)');
    ylabel('Spar spacing b (m)');
    title(['Skin Panel Buckling Stress (MPa) - ' labels{k}]);

    % 3D Surface Plot
    L_grid = reshape(L_vals_flat, 6, 6);
    b_grid = reshape(b_vals_flat, 6, 6);
    stress_grid = reshape(min_stress_vals, 6, 6);

    figure;
    surf(L_grid, b_grid, stress_grid);
    xlabel('Rib Spacing L (m)');
    ylabel('Spar Spacing b (m)');
    zlabel('Minimum Buckling Stress (MPa)');
    title(['3D Surface Plot - ' labels{k}]);
    shading interp;
    colormap jet;
    colorbar;

    % Line Plot: Stress vs b (for fixed L)
    figure; hold on;
    for i = 1:length(L_vals)
        idx = (L_vals_flat == L_vals(i));
        plot(b_vals_flat(idx), min_stress_vals(idx), '-o', ...
            'DisplayName', sprintf('L = %.2f m', L_vals(i)));
    end
    xlabel('Spar Spacing b (m)');
    ylabel('Minimum Buckling Stress (MPa)');
    title(['Min Stress vs Spar Spacing b - ' labels{k}]);
    legend('Location', 'northeastoutside');
    grid on;

    % Line Plot: Stress vs L (for fixed b)
    figure; hold on;
    for j = 1:length(b_vals)
        idx = (b_vals_flat == b_vals(j));
        plot(L_vals_flat(idx), min_stress_vals(idx), '-o', ...
            'DisplayName', sprintf('b = %.2f m', b_vals(j)));
    end
    xlabel('Rib Spacing L (m)');
    ylabel('Minimum Buckling Stress (MPa)');
    title(['Min Stress vs Rib Spacing L - ' labels{k}]);
    legend('Location', 'northeastoutside');
    grid on;
end

%% Summary Table Display for All Layups
for k = 1:3
    label = labels{k};
    layup = layups{k};
    Ncr_grid = Ncr_all{k};
    D = compute_D_matrix(layup, t_ply, E1, E2, G12, nu12);
    D11 = D(1,1); D12 = D(1,2); D22 = D(2,2); D66 = D(3,3);

    idx = 1;
    table_data = cell(36, 5);

    for i = 1:length(L_vals)
        for j = 1:length(b_vals)
            L = L_vals(i);
            b = b_vals(j);
            R = L / b;

            Ncr_vals = zeros(size(modes));
            for m = modes
                Ncr_vals(m) = (pi^2 / (m^2 * L^2)) * ...
                    (D11 * m^4 + 2*(D12 + 2*D66)*m^2 * R^2 + D22 * R^4);
            end

            [min_Ncr, best_m] = min(Ncr_vals);
            stress = min_Ncr / t_total / 1e6; % MPa

            table_data{idx,1} = L;
            table_data{idx,2} = b;
            table_data{idx,3} = stress;
            table_data{idx,4} = best_m;
            table_data{idx,5} = 'Skin Panel';
            idx = idx + 1;
        end
    end

    T_result = cell2table(table_data, ...
        'VariableNames', {'L_m', 'b_m', 'MinBucklingStress_MPa', 'Mode', 'FailedPart'});

    disp(['--- Results for Layup ' label ' ---']);
    disp(T_result);
end

%% Functions Used

function [w_tip, twist_tip] = compute_tip_response(EI, GJ, S, Kval, L, P)
% Computes tip deflection and twist using Rayleigh-Ritz method

Kmat = zeros(6,6);

% a-a terms (deflection)
Kmat(1,1) = 12*EI*L^3;
Kmat(1,2) = 36*EI*L^4;
Kmat(1,3) = 48*EI*L^5;
Kmat(2,2) = 144*EI*L^5 / 5;
Kmat(2,3) = 80*EI*L^6;
Kmat(3,3) = 400*EI*L^7 / 7;

% b-b terms (twist)
Kmat(4,4) = 12*S*L^3 + 9*GJ*L^5 / 5;
Kmat(4,5) = 36*S*L^4 + 4*GJ*L^6;
Kmat(4,6) = 48*S*L^5 + 30*GJ*L^7 / 7;
Kmat(5,5) = 144*S*L^5 / 5 + 16*GJ*L^7 / 7;
Kmat(5,6) = 80*S*L^6 + 5*GJ*L^8;
Kmat(6,6) = 400*S*L^7 / 7 + 25*GJ*L^9 / 9;

% a-b coupling terms (K)
Kab = zeros(3);
Kab(1,1) = 24*Kval*L^3;
Kab(1,2) = 36*Kval*L^4;
Kab(1,3) = 48*Kval*L^5;
Kab(2,1) = 36*Kval*L^4;
Kab(2,2) = 288*Kval*L^5 / 5;
Kab(2,3) = 80*Kval*L^6;
Kab(3,1) = 48*Kval*L^5;
Kab(3,2) = 80*Kval*L^6;
Kab(3,3) = 800*Kval*L^7 / 7;

% Fill in symmetric terms
Kmat(1:3, 4:6) = Kab;
Kmat(4:6, 1:3) = Kab';

% Load vector
f = zeros(6,1);
f(3) = P * L^5;

% Solve
q = Kmat \ f;
a = q(1:3);
b = q(4:6);

w_tip = a(1)*L^3 + a(2)*L^4 + a(3)*L^5;
twist_tip = b(1)*L^3 + b(2)*L^4 + b(3)*L^5;
end

function [a3, b3] = compute_tip_response_raw(EI, GJ, S, Kval, L, P)
% Returns a3 and b3 coefficients from Rayleigh-Ritz

Kmat = zeros(6,6);

% a-a terms
Kmat(1,1) = 12*EI*L^3;
Kmat(1,2) = 36*EI*L^4;
Kmat(1,3) = 48*EI*L^5;
Kmat(2,2) = 144*EI*L^5 / 5;
Kmat(2,3) = 80*EI*L^6;
Kmat(3,3) = 400*EI*L^7 / 7;

% b-b terms
Kmat(4,4) = 12*S*L^3 + 9*GJ*L^5 / 5;
Kmat(4,5) = 36*S*L^4 + 4*GJ*L^6;
Kmat(4,6) = 48*S*L^5 + 30*GJ*L^7 / 7;
Kmat(5,5) = 144*S*L^5 / 5 + 16*GJ*L^7 / 7;
Kmat(5,6) = 80*S*L^6 + 5*GJ*L^8;
Kmat(6,6) = 400*S*L^7 / 7 + 25*GJ*L^9 / 9;

% a-b coupling
Kab = zeros(3);
Kab(1,1) = 24*Kval*L^3;
Kab(1,2) = 36*Kval*L^4;
Kab(1,3) = 48*Kval*L^5;
Kab(2,1) = 36*Kval*L^4;
Kab(2,2) = 288*Kval*L^5 / 5;
Kab(2,3) = 80*Kval*L^6;
Kab(3,1) = 48*Kval*L^5;
Kab(3,2) = 80*Kval*L^6;
Kab(3,3) = 800*Kval*L^7 / 7;

Kmat(1:3, 4:6) = Kab;
Kmat(4:6, 1:3) = Kab';

f = zeros(6,1);
f(3) = P * L^5;

q = Kmat \ f;
a = q(1:3);
b = q(4:6);

a3 = a(3);
b3 = b(3);
end

function Dmat = compute_D_matrix(ply_angles_deg, t_ply, E1, E2, G12, nu12)
% Computes the bending stiffness D matrix for a laminate

N = length(ply_angles_deg);
total_thickness = N * t_ply;
nu21 = (nu12 * E2) / E1;

S = [1/E1, -nu12/E1, 0;
    -nu12/E1, 1/E2, 0;
    0, 0, 1/G12];

Dmat = zeros(3,3);
z = linspace(-total_thickness/2, total_thickness/2, N+1);

for k = 1:N
    theta = deg2rad(ply_angles_deg(k));
    c = cos(theta);
    s = sin(theta);

    T_sigma = [c^2, s^2, 2*c*s;
               s^2, c^2, -2*c*s;
              -c*s, c*s, c^2 - s^2];

    T_epsilon = [c^2, s^2, c*s;
                 s^2, c^2, -c*s;
                -2*c*s, 2*c*s, c^2 - s^2];

    Q = inv(T_sigma) * inv(S) * inv(T_epsilon);

    z_bot = z(k);
    z_top = z(k+1);
    delta_z3 = (z_top^3 - z_bot^3) / 3;

    Dmat = Dmat + Q * delta_z3;
end
end

function Q = get_Q_matrix(E1, E2, G12, nu12)
% Returns the reduced stiffness matrix Q

nu21 = (nu12 * E2) / E1;
denom = 1 - nu12 * nu21;

Q11 = E1 / denom;
Q12 = nu12 * E2 / denom;
Q22 = E2 / denom;
Q66 = G12;

Q = [Q11, Q12, 0;
     Q12, Q22, 0;
     0, 0, Q66];
end

function Qbar = transform_Q(Q, theta_deg)
% Returns the transformed stiffness matrix Qbar for a given angle

theta = deg2rad(theta_deg);
c = cos(theta);
s = sin(theta);

T_sigma = [c^2, s^2, 2*c*s;
           s^2, c^2, -2*c*s;
          -c*s, c*s, c^2 - s^2];

T_epsilon = [c^2, s^2, c*s;
             s^2, c^2, -c*s;
            -2*c*s, 2*c*s, c^2 - s^2];

Qbar = T_sigma \ (Q * T_epsilon);
end

%% Helper Integrals for Rayleigh Quotient

function val = integral_phi_dd_squared(L)
% Computes integral of (d^2/dx^2(x^3))^2 from 0 to L
val = 36 * integral(@(x) x.^2, 0, L);  % 36 * L^3 / 3 = 12 * L^3
end

function val = integral_phi_squared(L)
% Computes integral of (x^3)^2 from 0 to L
val = L^7 / 7;
end

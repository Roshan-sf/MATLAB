%% Problem 1 - Parts A through E
% Geometry and load
L = 1.0; % Length in meters
b = 0.1; % Width in meters
t_total = 0.025; % Total thickness in meters
P = 1000; % Tip load in Newtons
% Ply data
Nplies = 5;
t_ply = t_total / Nplies;
% Material properties for Carbon/Epoxy
E1 = 140e9; % Pa
E2 = 10e9; % Pa
G12 = 7e9; % Pa
nu12 = 0.3;
% Part A
% Compute tip deflection and twist for five symmetric layups using the
% Rayleigh-Ritz method.
% Defining layups (bottom ply is positive)
layup_angles = {
[75 -75 75 -75 75],
[60 -60 60 -60 60],
[45 -45 45 -45 45],
[30 -30 30 -30 30],
[10 -10 10 -10 10]
};
labels = {'[±75°]', '[±60°]', '[±45°]', '[±30°]', '[±10°]'};
% Initializing result arrays
w_tips = zeros(1,5);
twist_tips = zeros(1,5);
disp('----Part A----')
for i = 1:5
theta = layup_angles{i};
D = compute_D_matrix(theta, t_ply, E1, E2, G12, nu12);
% Extracting stiffness constants
D11 = D(1,1); D12 = D(1,2); D16 = D(1,3);
D66 = D(3,3);
EI = D11 * b;
GJ = 4 * D66 * b;
K = 2 * D16 * b;
S = (1/12) * D11 * b^3;
[w_tip, twist_tip] = compute_tip_response(EI, GJ, S, K, L, P);
w_tips(i) = w_tip;
twist_tips(i) = twist_tip;
% Converting twist to degrees for display
twist_deg = twist_tip * (180/pi);
disp([labels{i} ': w_tip = ' num2str(w_tip, '%.4e') ' m, twist = '
num2str(twist_deg, '%.4e') ' deg']);
end
% Part B
% Design and analyze a symmetric layup that results in zero tip twist
% under loading.
% Custom layup testing: designed to reduce twist
% layup = [60 -30 0 -30 60]; (Doesnt work)
% layup = [30 -60 0 -60 30]; (Doesnt work)
% layup = [10 -45 0 -45 10]; (Doesnt work)
layup = [0 0 0 0 0]; % Works!
% Computing ABD matrix and extract D matrix terms
D = compute_D_matrix(layup, t_ply, E1, E2, G12, nu12);
D11 = D(1,1); D12 = D(1,2); D16 = D(1,3);
D66 = D(3,3);
% Computing stiffness terms
EI = D11 * b;
GJ = 4 * D66 * b;
K = 2 * D16 * b;
S = (1/12) * D11 * b^3;
% Computing tip response
[w_tip, twist_tip] = compute_tip_response(EI, GJ, S, K, L, P);
twist_deg = twist_tip * (180/pi);
% Displaying results
disp('----Part B----')
disp(['Custom Layup: ' num2str(layup)]);
disp(['w_tip = ' num2str(w_tip, '%.4e') ' m']);
disp(['twist = ' num2str(twist_deg, '%.4e') ' deg']);
% Displaying a recommendation if twist is nearly zero for testing
if abs(twist_deg) < 1e-3
disp('Twist is effectively zero — this layup is a strong candidate.');
else
disp('Twist is not zero — consider modifying layup or testing others.');
end
% Part C
% Evaluate five new symmetric layups to identify configurations that
% minimize tip deflection.
% Defining 5 new symmetric layups
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
% Initializing result arrays
w_tips = zeros(1,5);
twist_tips = zeros(1,5);
disp('----Part C----')
for i = 1:length(layups)
layup = layups{i};
label = labels{i};
% Computing D matrix
D = compute_D_matrix(layup, t_ply, E1, E2, G12, nu12);
D11 = D(1,1); D12 = D(1,2); D16 = D(1,3);
D66 = D(3,3);
% Computing stiffness terms
EI = D11 * b;
GJ = 4 * D66 * b;
K = 2 * D16 * b;
S = (1/12) * D11 * b^3;
% Computing response
[w_tip, twist_tip] = compute_tip_response(EI, GJ, S, K, L, P);
twist_deg = twist_tip * (180/pi);
% Storing results
w_tips(i) = w_tip;
twist_tips(i) = twist_deg;
% Displaying results
disp(['Layup ' label ':']);
disp(['w_tip = ' num2str(w_tip, '%.4e') ' m']);
disp(['twist = ' num2str(twist_deg, '%.4e') ' deg']);
disp('');
end
% Part D
% Identify the layup with the highest twist and determine the tip load at
% which ply failure occurs based on the Tsai-Hill failure index.
% More material properties (Carbon/Epoxy)
S1 = 1500e6;
S2 = 40e6;
S12 = 70e6; % Strengths in Pa
% Layup with max twist
layup = [0 -60 0 -60 0];
label = '[0 -60 0 -60 0]';
% Computing D matrix and curvature coefficients
D = compute_D_matrix(layup, t_ply, E1, E2, G12, nu12);
D11 = D(1,1); D12 = D(1,2); D16 = D(1,3);
D66 = D(3,3);
EI = D11 * b;
GJ = 4 * D66 * b;
K = 2 * D16 * b;
S = (1/12) * D11 * b^3;
% Using curvature to compute strain per ply for increasing load
load_vals = linspace(50, 5000, 500); % Try tip loads from 50 to 5000 N
% Getting z positions of ply midplanes (symmetric stack)
z_mids = linspace(-t_total/2 + t_ply/2, t_total/2 - t_ply/2, Nplies);
failure_load = NaN;
failure_ply = NaN;
failure_index = NaN;
for P = load_vals
[a3, b3] = compute_tip_response_raw(EI, GJ, S, K, L, P); % use version
returning a3, b3
kappa = [a3; b3; 2*a3];
for i = 1:Nplies
theta = layup(i);
z = z_mids(i);
eps = z * kappa; % midplane strains
% Qbar for each ply
Q = get_Q_matrix(E1, E2, G12, nu12);
Qbar = transform_Q(Q, theta);
sigma = Qbar * eps;
sigma1 = sigma(1); sigma2 = sigma(2); tau12 = sigma(3);
% Tsai-Hill failure index
FI = (sigma1/S1)^2 - (sigma1*sigma2)/(S1^2) + (sigma2/S2)^2 + (tau12/S12)^2;
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
% Reporting results
disp('----Part D----')
disp(['Layup with max twist: ' label]);
if ~isnan(failure_load)
disp(['Failure occurs at P = ' num2str(failure_load, '%.2f') ' N']);
disp(['Critical ply index: ' num2str(failure_ply)]);
disp(['Failure index: ' num2str(failure_index, '%.3f')]);
else
disp('No failure detected up to maximum load.');
end
% Part E
% Calculate and compare the lowest natural frequencies of all 11 layups
% using the Rayleigh quotient method.
% Mass property
rho = 1615; % density of carbon/epoxy in kg/m^3
% Layups (from parts A–C)
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
labels = {'[±75°]', '[±60°]', '[±45°]',...
'[±30°]', '[±10°]', '[0°]', '[0/15]',...
'[±10/0]', '[0/30]', '[0/-45]', '[0/-60]'};
frequencies_Hz = zeros(1, length(layups));
disp('----Part E----')
for i = 1:length(layups)
theta = layups{i};
D = compute_D_matrix(theta, t_ply, E1, E2, G12, nu12);
D11 = D(1,1); D66 = D(3,3);
EI = D11 * b;
GJ = 4 * D66 * b;
% Approximating mass
m = rho * b * t_total * L; % total mass
% Using Rayleigh Quotient for fundamental frequency (bending only)
% Mode shape: phi(x) = x^3 (consistent with h(x))
% w = sqrt((EI * integral(phi''^2)) / (m * integral(phi^2)))
num = EI * integral_phi_dd_squared(L);
den = m * integral_phi_squared(L);
omega_rad = sqrt(num / den); % rad/s
freq_Hz = omega_rad / (2*pi); % Hz
frequencies_Hz(i) = freq_Hz;
disp([labels{i} ': f = ' num2str(freq_Hz, '%.2f') ' Hz']);
end
%% Problem 2
% Buckling of composite skin panel using laminate theory
% Material properties for graphite/epoxy composite
E1 = 181e9; % Pa
E2 = 10.3e9; % Pa
G12 = 7.17e9; % Pa
nu12 = 0.28;
% Geometry
t_total = 0.005; % Total laminate thickness [m] (5 mm)
Nplies = 5;
t_ply = t_total / Nplies;
% Grid (same as Assignment 4 Q1)
L_vals = linspace(0.1, 0.5, 6); % Rib spacing [m]
b_vals = linspace(0.1, 0.3, 6); % Spar spacing [m]
% Mode numbers to check
modes = 1:6;
% Defining symmetric layups
layups = {
[0 0 0 0 0], % [0]_s
[0 90 90 90 0], % [0/90]_s
[90 0 0 0 90] % [90/0]_s
};
labels = {'[0]_s', '[0/90]_s', '[90/0]_s'};
% Preallocating results
Ncr_all = cell(1,3);
for k = 1:3
layup = layups{k};
Ncr_grid = zeros(length(L_vals), length(b_vals));
min_stress_vals = zeros(36,1);
L_vals_flat = zeros(36,1);
b_vals_flat = zeros(36,1);
% Computing D matrix for this layup
D = compute_D_matrix(layup, t_ply, E1, E2, G12, nu12);
D11 = D(1,1); D12 = D(1,2); D22 = D(2,2); D66 = D(3,3);
idx_flat = 1;
for i = 1:length(L_vals)
for j = 1:length(b_vals)
L = L_vals(i);
b = b_vals(j);
R = L/b;
% Looping over m and compute N_cr for each
Ncr_vals = zeros(size(modes));
for m = modes
Ncr_vals(m) = (pi^2 / (m^2 * L^2)) * ...
(D11 * m^4 + 2*(D12 + 2*D66)*m^2 * R^2 + D22 * R^4);
end
% Taking minimum across m
Ncr_grid(i,j) = min(Ncr_vals);
min_stress_vals(idx_flat) = Ncr_grid(i,j) / t_total / 1e6; % Convert to
MPa
L_vals_flat(idx_flat) = L;
b_vals_flat(idx_flat) = b;
idx_flat = idx_flat + 1;
end
end
Ncr_all{k} = Ncr_grid;
% Skin Panel Contour Plot
figure;
contourf(L_vals, b_vals, Ncr_grid' / t_total / 1e6, 20, 'LineColor', 'none');
colorbar;
xlabel('Rib spacing L (m)'); ylabel('Spar spacing b (m)');
title(['Skin Panel Buckling Stress (MPa) - ' labels{k}]);
% Additional Plots for Layup labels{k}
% Reshaping grids
L_grid = reshape(L_vals_flat, 6, 6);
b_grid = reshape(b_vals_flat, 6, 6);
stress_grid = reshape(min_stress_vals, 6, 6); % MPa
% 3D Surface Plot
figure;
surf(L_grid, b_grid, stress_grid);
xlabel('Rib Spacing L (m)');
ylabel('Spar Spacing b (m)');
zlabel('Minimum Buckling Stress (MPa)');
title(['3D Surface Plot of Min Stress - ' labels{k}]);
shading interp;
colormap jet;
colorbar;
% Stress vs b for fixed L
figure;
hold on;
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
% Stress vs L for fixed b
figure;
hold on;
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
% Displaying Summary Tables for All Layups
for k = 1:3
label = labels{k};
Ncr_grid = Ncr_all{k};
idx = 1;
table_data = cell(36, 5);
% Recomputing D matrix for this layup
layup = layups{k};
D = compute_D_matrix(layup, t_ply, E1, E2, G12, nu12);
D11 = D(1,1); D12 = D(1,2); D22 = D(2,2); D66 = D(3,3);
for i = 1:length(L_vals)
for j = 1:length(b_vals)
L = L_vals(i);
b = b_vals(j);
R = L / b;
% Looping over m and compute N_cr for each
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
table_data{idx,4} = best_m; % Mode number giving min stress
table_data{idx,5} = 'Skin Panel';
idx = idx + 1;
end
end
T_result = cell2table(table_data, ...
'VariableNames', {'L_m', 'b_m', 'MinBucklingStress_MPa', 'Mode',
'FailedPart'});
disp(['--- Results for Layup ' label ' ---']);
disp(T_result);
end
%% Functions used
function [w_tip, twist_tip] = compute_tip_response(EI, GJ, S, Kval, L, P)
% Computes tip deflection and twist for a composite cantilever beam using Rayleigh-
Ritz
% Symbolic stiffness matrix entries (from integral of energy)
Kmat = zeros(6,6);
% a-a terms (deflection)
Kmat(1,1) = 12*EI*L^3;
Kmat(1,2) = 36*EI*L^4;
Kmat(1,3) = 48*EI*L^5;
Kmat(2,2) = 144*EI*L^5/5;
Kmat(2,3) = 80*EI*L^6;
Kmat(3,3) = 400*EI*L^7/7;
% b-b terms (twist)
Kmat(4,4) = 12*S*L^3 + 9*GJ*L^5/5;
Kmat(4,5) = 36*S*L^4 + 4*GJ*L^6;
Kmat(4,6) = 48*S*L^5 + 30*GJ*L^7/7;
Kmat(5,5) = 144*S*L^5/5 + 16*GJ*L^7/7;
Kmat(5,6) = 80*S*L^6 + 5*GJ*L^8;
Kmat(6,6) = 400*S*L^7/7 + 25*GJ*L^9/9;
% a-b coupling terms (K)
Kab = zeros(3);
Kab(1,1) = 24*Kval*L^3;
Kab(1,2) = 36*Kval*L^4;
Kab(1,3) = 48*Kval*L^5;
Kab(2,1) = 36*Kval*L^4;
Kab(2,2) = 288*Kval*L^5/5;
Kab(2,3) = 80*Kval*L^6;
Kab(3,1) = 48*Kval*L^5;
Kab(3,2) = 80*Kval*L^6;
Kab(3,3) = 800*Kval*L^7/7;
% Filling in symmetric matrix
Kmat(1:3, 4:6) = Kab;
Kmat(4:6, 1:3) = Kab';
% Force vector (only last deflection mode has non-zero contribution at tip)
f = zeros(6,1);
f(3) = P * L^5; % h(x) = a3 x^5 => w_tip = a3 L^5
% Solving the system
q = Kmat \ f;
% Extracting tip deflection and twist
a = q(1:3);
b = q(4:6);
w_tip = a(1)*L^3 + a(2)*L^4 + a(3)*L^5;
twist_tip = b(1)*L^3 + b(2)*L^4 + b(3)*L^5;
end
function Dmat = compute_D_matrix(ply_angles_deg, t_ply, E1, E2, G12, nu12)
% Number of plies
N = length(ply_angles_deg);
total_thickness = N * t_ply;
% Minor Poisson's ratio
nu21 = (nu12 * E2) / E1;
% Compliance matrix S for a single ply
S = [1/E1, -nu12/E1, 0;
-nu12/E1, 1/E2, 0;
0, 0, 1/G12];
% Initialize D matrix
Dmat = zeros(3,3);
% z-coordinates (bottom to top)
z = linspace(-total_thickness/2, total_thickness/2, N+1);
for k = 1:N
theta = ply_angles_deg(k) * pi / 180; % convert to radians
c = cos(theta);
s = sin(theta);
% Transformation matrices
T_sigma = [c^2, s^2, 2*c*s;
s^2, c^2, -2*c*s;
-c*s, c*s, c^2 - s^2];
T_epsilon = [c^2, s^2, c*s;
s^2, c^2, -c*s;
-2*c*s, 2*c*s, c^2 - s^2];
Q = inv(T_sigma) * inv(S) * inv(T_epsilon); % Transformed Q-bar
z_top = z(k+1);
z_bot = z(k);
delta_z3 = (z_top^3 - z_bot^3) / 3;
% Accumulate D matrix
Dmat = Dmat + Q * delta_z3;
end
end
function [a3, b3] = compute_tip_response_raw(EI, GJ, S, Kval, L, P)
% Returns a3 and b3 from Rayleigh-Ritz solution for deflection and twist
% Stiffness matrix (same as original)
Kmat = zeros(6,6);
% a-a terms (deflection)
Kmat(1,1) = 12*EI*L^3;
Kmat(1,2) = 36*EI*L^4;
Kmat(1,3) = 48*EI*L^5;
Kmat(2,2) = 144*EI*L^5/5;
Kmat(2,3) = 80*EI*L^6;
Kmat(3,3) = 400*EI*L^7/7;
% b-b terms (twist)
Kmat(4,4) = 12*S*L^3 + 9*GJ*L^5/5;
Kmat(4,5) = 36*S*L^4 + 4*GJ*L^6;
Kmat(4,6) = 48*S*L^5 + 30*GJ*L^7/7;
Kmat(5,5) = 144*S*L^5/5 + 16*GJ*L^7/7;
Kmat(5,6) = 80*S*L^6 + 5*GJ*L^8;
Kmat(6,6) = 400*S*L^7/7 + 25*GJ*L^9/9;
% a-b coupling terms (K)
Kab = zeros(3);
Kab(1,1) = 24*Kval*L^3;
Kab(1,2) = 36*Kval*L^4;
Kab(1,3) = 48*Kval*L^5;
Kab(2,1) = 36*Kval*L^4;
Kab(2,2) = 288*Kval*L^5/5;
Kab(2,3) = 80*Kval*L^6;
Kab(3,1) = 48*Kval*L^5;
Kab(3,2) = 80*Kval*L^6;
Kab(3,3) = 800*Kval*L^7/7;
Kmat(1:3, 4:6) = Kab;
Kmat(4:6, 1:3) = Kab';
% Load vector
f = zeros(6,1);
f(3) = P * L^5;
% Solve
q = Kmat \ f;
a = q(1:3);
b = q(4:6);
a3 = a(3);
b3 = b(3);
end
function Q = get_Q_matrix(E1, E2, G12, nu12)
% Minor Poisson's ratio
nu21 = (nu12 * E2) / E1;
% Denominator for stiffness terms
denom = 1 - nu12 * nu21;
% Reduced stiffness matrix Q
Q11 = E1 / denom;
Q12 = nu12 * E2 / denom;
Q22 = E2 / denom;
Q66 = G12;
Q = [Q11, Q12, 0;
Q12, Q22, 0;
0, 0, Q66];
end
function Qbar = transform_Q(Q, theta_deg)
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
% Helper integrals for Rayleigh Quotient
function val = integral_phi_dd_squared(L)
% Integral of (d^2/dx^2 (x^3))^2 from 0 to L => (6x)^2 = 36x^2
val = 36 * integral(@(x) x.^2, 0, L); % = 36*L^3/3 = 12L^3
end
function val = integral_phi_squared(L)
% Integral of (x^3)^2 = x^6 from 0 to L => L^7 / 7
val = L^7 / 7;
end
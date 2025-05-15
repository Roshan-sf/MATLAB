%% Roshan Jaiswal-Ferri, Alessandro Tedeschi, Stefan Rosu
%Section - 02
%Aero 431 HW2: 4/28/25

%% Workspace Prep

%warning off
format long     %Allows for more accurate decimals
close all;      %Clears all
clear all;      %Clears Workspace
clc;            %Clears Command 

%% Rayleigh-Ritz 

% Given constants
E = 70e9; % Young's modulus [Pa]
nu = 0.3; % Poisson's ratio
rho = 2710; % Density [kg/m^3]
h = 0.01; % Thickness [m]
L = 1; % Plate side length [m]
a_val = 0; % Assume a = 0 for second ansatz in part (a)

% Flexural rigidity
D = E * h^3 / (12 * (1 - nu^2));

% First mode shape
% w1 = (x^2 - L^2/4)^2 * (y^2 - L^2/4)

% Symbolic integration values (computed earlier)
num1 = 72 * E * h^2;
den1 = L^4 * rho * (1 - nu^2);
omega1_sq = num1 / den1;
omega1_rad = sqrt(omega1_sq); % in rad/s
freq1 = omega1_rad / (2*pi); % in Hz

% Second Ansatz with a = 0
% w2 = (x^2 - L^2/4)^2 * (y^2 - L^2/4)^2 * (1 + a*(x^2 + y^2))

num2 = 1584 * E * h^2 * 429;
den2 = L^4 * rho * (6292 * (1 - nu^2));
omega2_sq = num2 / den2;
omega2_rad = sqrt(omega2_sq);
freq2 = omega2_rad / (2*pi);

% Displaying Results
disp('Estimated Fundamental Frequencies:');
disp(['First mode shape (simple): ', num2str(freq1), ' Hz']);
disp(['Second mode shape (refined, a=0): ', num2str(freq2), ' Hz']);

%% Problem 2 - Part b

% Constants
E = 70e9; % Young's modulus [Pa]
nu = 0.3; % Poisson's ratio
rho = 2710; % Density [kg/m^3]
h = 0.01; % Plate thickness [m]
L = 1.0; % Plate side length [m]

% Flexural rigidity and mass per area
D = E * h^3 / (12 * (1 - nu^2));
rhoh = rho * h;

% Mass matrix M (4x4)
M = [L^10/900, 0, 0, L^12/12600;
0, L^12/25200, 0, 0;
0, 0, L^12/25200, 0;
L^12/12600, 0, 0, L^14/105840];

% Stiffness matrix K (from symbolic integration)
K = [ ...
    0.0013990, 0, 0, 7.1737e-04;
    0, 3.7919e-04, 0, 0;
    0, 0, 3.7919e-04, 0;
    7.1737e-04, 0, 0, 2.3851e-04];

% Scale by D and rho*h
K_scaled = D * K;
M_scaled = rhoh * M;

% Solving generalized eigenvalue problem
[mode_shapes, omega_sq] = eig(K_scaled, M_scaled);
omega_sq_vals = diag(omega_sq);
omega_vals = sqrt(omega_sq_vals);
frequencies = omega_vals / (2*pi); % Convert to Hz

% Sorting and displaying first 4 frequencies
[frequencies_sorted, idx] = sort(real(frequencies));
disp('Estimated Natural Frequencies (Hz) from Generalized Eigenvalue Problem:');
for i = 1:4
    disp(['Mode ', num2str(i), ': ', num2str(frequencies_sorted(i), '%.4f'),...
        'Hz']);
end

% Plotting Mode Shapes using the original ansatz function
% w(x,y) = weight(x,y) * (a + bx + cy + d(x^2 + y^2))

phi1 = @(x,y) 1;
phi2 = @(x,y) x;
phi3 = @(x,y) y;
phi4 = @(x,y) x.^2 + y.^2;

wgt = @(x,y) (x.^2 - L^2/4).^2 .* (y.^2 - L^2/4).^2;

[X, Y] = meshgrid(linspace(-L/2, L/2, 100), linspace(-L/2, L/2, 100));

for mode_num = 1:4
    coeffs = mode_shapes(:, idx(mode_num));
    W = wgt(X, Y) .* ( ...
    coeffs(1)*phi1(X,Y) + ...
    coeffs(2)*phi2(X,Y) + ...
    coeffs(3)*phi3(X,Y) + ...
    coeffs(4)*phi4(X,Y));
    figure;
    surf(X, Y, W, 'EdgeColor', 'none');
    colormap jet;
    colorbar;
    title(['Mode ', num2str(mode_num), ' Shape (Hz: ',...
    num2str(frequencies_sorted(mode_num), '%.2f'), ')']);
    xlabel('x [m]');
    ylabel('y [m]');
    zlabel('w(x,y)');
    view(30, 30);
    axis tight;
end
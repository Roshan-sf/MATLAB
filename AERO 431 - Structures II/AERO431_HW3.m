%% Roshan Jaiswal-Ferri, Alessandro Tedeschi, Stefan Rosu
%Section - 02
%Aero 431 HW3: 5/12/25

%% Workspace Prep

%warning off
format long     %Allows for more accurate decimals
close all;      %Clears all
clear all;      %Clears Workspace
clc;            %Clears Command 

%% Rayleigh-Ritz 

% Material and plate properties
E = 70e9;         % Youngs Modulus Pa
nu = 0.3;         % Poisson's Ratio
rho = 2710;       % kg/m^3
h = 0.01;         % thickness m
L = 1;          % side length m

% rigidity
D = E*h^3 / (12*(1-nu^2));

%% Part 1 — Rayleigh Quotient Method

% symbolic Rayleigh quotient
omega1_sq = (72*E*h^2) / (L^4*rho*(1 - nu^2));
freq1 = sqrt(omega1_sq) / (2*pi);  % Hz

omega2_sq = (1584*E*h^2*429) / (L^4*rho*6292*(1-nu^2));
freq2 = sqrt(omega2_sq) / (2*pi);  % Hz

fprintf('--- Fundamental Frequency Estimates (Rayleigh Method) ---\n');
fprintf('1st Ansatz (simple):      %.4f Hz\n', freq1);
fprintf('2nd Ansatz (refined, a=0): %.4f Hz\n\n', freq2);

%% Part 2 — Generalized Eigenvalues

% Mass / area
rhoh = rho*h;

% mass and stiffness matrix
M = [L^10/900,     0,               0,           L^12/12600;
     0,           L^12/25200,      0,           0;
     0,           0,               L^12/25200,  0;
     L^12/12600,  0,               0,           L^14/105840];

K = [0.0013990,     0,          0,          7.1737e-04;
     0,             3.7919e-04, 0,          0;
     0,             0,          3.7919e-04, 0;
     7.1737e-04,    0,          0,          2.3851e-04];

K_scaled = D * K;
M_scaled = rhoh * M;

% solve generalized eigenvalue problem
[mode_shapes, omega_sq_vals] = eig(K_scaled, M_scaled);
omega_vals = sqrt(diag(omega_sq_vals));
frequencies = omega_vals / (2 * pi);

% Sort frequencies
[frequencies_sorted, sort_idx] = sort(real(frequencies));

fprintf('--- Frequencies from Generalized Eigenvalue Problem ---\n');
for i = 1:4
    fprintf('Mode %d: %.4f Hz\n', i, frequencies_sorted(i));
end

% Comments:
% The simple mode found a frequency about 16hz lower than the refined
% mode. The inclusion of more terms for displacement shows that the refined
% mode is creating a better aproximation of the actual frequency. 

%% Plotting

phi = {
    @(x, y) 1;
    @(x, y) x;
    @(x, y) y;
    @(x, y) x.^2 + y.^2
};

% weight function
wgt = @(x, y) (x.^2 - L^2/4).^2 .* (y.^2 - L^2/4).^2;

x_vals = linspace(-L/2, L/2, 120);
y_vals = linspace(-L/2, L/2, 80);
[Y, X] = ndgrid(y_vals, x_vals);

for i = 1:4
    coeffs = mode_shapes(:, sort_idx(i));
    
    W = wgt(X, Y);
    for j = 1:4
        W = W + coeffs(j) * wgt(X, Y) .* phi{j}(X, Y);
    end
    
    % Plot
    figure;
    surf(X, Y, W, 'EdgeColor', 'none');
    colormap turbo;
    colorbar;
    title(sprintf('Mode %d Shape (%.2f Hz)', i, frequencies_sorted(i)));
    xlabel('x [m]');
    ylabel('y [m]');
    zlabel('w(x,y)');
    view(40, 25);
    axis tight;
end

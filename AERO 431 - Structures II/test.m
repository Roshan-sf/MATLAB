%% Roshan Jaiswal-Ferri, Alessandro Tedeschi, Stefan Rosu
%Section - 02
%Aero 431 HW2: 4/28/25

%% Rayleigh-Ritz Method using Quadratic, Cubic, Quartic, and Quintic Trial Functions

E = 70e9;    % Pa
A = 5e-6;    % m^2
L = 1;       % m
P = 100e3;   % N/m

syms x real

% Exact solution
u_exact_fn = @(xq) (P*L^2/(4*E*A)) * (xq/L - (1/5)*(xq/L).^5);

degrees = [2, 3, 4, 5];
colors = {'b--', 'r-.', 'g:', 'm-.'};
labels = {'Quadratic', 'Cubic', 'Quartic', 'Quintic'};

% Plot setup
x_plot = linspace(0, L, 500);
figure;
hold on;
plot(x_plot, u_exact_fn(x_plot), 'k-', 'LineWidth', 2, 'DisplayName', 'Exact Solution');

for idx = 1:length(degrees)
    n = degrees(idx);
    
    % Define trial function manually
    a = sym('a', [n,1], 'real'); % n coefficients: a1, a2, ..., an
    
    u_trial = 0;
    for i = 1:n
        u_trial = u_trial + a(i)*x^(i); % build a1*x^1 + a2*x^2 + a3*x^3 + ...
    end
    % Make sure it satisfies u(0) = 0: already done (no constant term)
    
    % Compute energy functional
    u_trial_dx = diff(u_trial, x);
    
    % Potential energy functional E(u)
    b = -P*(x/L)^3;  % <<<< FLIPPED SIGN here
    Energy = (1/2)*E*A*int(u_trial_dx^2, x, 0, L) - int(b*u_trial, x, 0, L);
    
    % Derive system of equations
    eqns = gradient(Energy, a);
    eqns = simplify(eqns);
    [M, f] = equationsToMatrix(eqns, a);
    
    M = double(M);
    f = double(-f);
    
    % Solve for coefficients
    a_sol = M\f;
    
    % Substitute coefficients back into u_trial
    u_ritz_expr = subs(u_trial, a, a_sol);
    u_ritz_fn = matlabFunction(u_ritz_expr, 'Vars', x);
    
    % Plot
    plot(x_plot, u_ritz_fn(x_plot), colors{idx}, 'LineWidth', 2, 'DisplayName', labels{idx});
end

xlabel('x (m)');
ylabel('Displacement u(x) (m)');
title('Rayleigh-Ritz Approximation (Quadratic to Quintic)');
legend('Location', 'northwest');
grid on;
hold off;


%% FEM Solution for 3, 4, and 5 elements using piecewise linear hat functions

% Problem parameters
E = 70e9;    % Young's Modulus [Pa]
A = 5e-6;    % Cross-sectional Area [m^2]
L = 1;       % Length of rod [m]
P = 100e3;   % Load coefficient [N/m]

% Exact solution
u_exact_fn = @(xq) (P*L^2/(4*E*A)) * (xq/L - (1/5)*(xq/L).^5);

% Solve and plot for 3, 4, and 5 elements
n_list = [3, 4, 5];
colors = {'b--', 'r-.', 'g:'};
x_plot = linspace(0, L, 500);

figure;
hold on;
plot(x_plot, u_exact_fn(x_plot), 'k-', 'LineWidth', 2, 'DisplayName', 'Exact');

for idx = 1:length(n_list)
    n_elements = n_list(idx);
    [x_nodes, U] = FEM_solver(n_elements, E, A, L, P);
    
    % Interpolated FEM solution
    u_fem_fn = @(xq) fem_interpolate(xq, x_nodes, U);
    
    plot(x_plot, arrayfun(u_fem_fn, x_plot), colors{idx}, 'LineWidth', 2, ...
         'DisplayName', sprintf('FEM (%d elements)', n_elements));
end

xlabel('x (m)');
ylabel('Displacement u(x) (m)');
title('Finite Element Solution vs Exact Solution');
legend('Location', 'northwest');
grid on;
hold off;

%% FEM Solver Function
function [x_nodes, U_full] = FEM_solver(n_elements, E, A, L, P)
    syms x real
    
    % Mesh
    x_nodes = linspace(0, L, n_elements+1);
    h = x_nodes(2) - x_nodes(1);
    
    % Stiffness matrix K and load vector F (full size, including all nodes)
    K_full = zeros(n_elements+1);
    F_full = zeros(n_elements+1,1);
    
    % Loop over elements
    for e = 1:n_elements
        x1 = x_nodes(e);
        x2 = x_nodes(e+1);
        
        % Local stiffness matrix for element
        Ke = (E*A/h) * [1 -1; -1 1];
        
        % Local load vector
        b = @(xq) P*(xq/L).^3;  % <-- FIXED HERE (. instead of ^)
        phi1 = @(xq) (x2 - xq)/h;
        phi2 = @(xq) (xq - x1)/h;
        
        Fe = zeros(2,1);
        Fe(1) = integral(@(xq) b(xq).*phi1(xq), x1, x2);
        Fe(2) = integral(@(xq) b(xq).*phi2(xq), x1, x2);
        
        % Assembly
        idx = [e, e+1];
        K_full(idx, idx) = K_full(idx, idx) + Ke;
        F_full(idx) = F_full(idx) + Fe;
    end
    
    % Apply boundary condition u(0) = 0
    K_reduced = K_full(2:end, 2:end);
    F_reduced = F_full(2:end);
    
    % Solve for unknown displacements
    U_reduced = K_reduced \ F_reduced;
    
    % Insert known boundary displacement u(0) = 0
    U_full = [0; U_reduced];
end

%% FEM Interpolation Function
function u_val = fem_interpolate(xq, x_nodes, U)
    n_elements = length(x_nodes) - 1;
    h = x_nodes(2) - x_nodes(1);
    u_val = zeros(size(xq));
    
    for e = 1:n_elements
        x1 = x_nodes(e);
        x2 = x_nodes(e+1);
        idx = (xq >= x1) & (xq <= x2);
        if any(idx)
            xi = xq(idx);
            phi1 = (x2 - xi)/h;
            phi2 = (xi - x1)/h;
            u_val(idx) = U(e)*phi1 + U(e+1)*phi2;
        end
    end
end

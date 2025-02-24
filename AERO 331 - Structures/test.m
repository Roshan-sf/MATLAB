clc; clear; close all;

% Given dimensions
h = 0.1;        % Height of the section (m)
b = 0.2;        % Width of the section (m)
t = 0.005;      % Thickness of the walls (m)
L_range = linspace(0.025, 0.175, 100); % L values from 0.025m to 0.175m

% Compute torsion constant for single-cell section
A1 = h * b; % Enclosed area
perimeter1 = 2*(h + b); % Perimeter
J1 = (4 * A1^2) / (perimeter1 / t); % Bredt's formula

% Compute torsion constant for two-cell section
J2 = zeros(size(L_range));
for i = 1:length(L_range)
    L = L_range(i);
    A2_left = h * L;
    A2_right = h * (b - L);
    
    % Perimeters
    s1 = 2 * (L + h);
    s2 = 2 * (b - L + h);
    s3 = h; % The new web
    
    % Torsion constant for the two-cell section (compatibility method)
    J2(i) = (4 * A2_left^2) / (s1 / t) + (4 * A2_right^2) / (s2 / t) + (4 * (A2_left * A2_right) / s3) / t;
end

% Compute the ratio of torsion constants
ratio = J2 / J1;

% Plot the results
figure;
plot(L_range, ratio, 'b', 'LineWidth', 2);
grid on;
xlabel('L (m)');
ylabel('J_{two-cell} / J_{single-cell}');
title('Ratio of Two-Cell to Single-Cell Torsion Constant');

% Davis Sok
% Aero 300 LAB5 : 2024.05.08

%close all;              % clears all
clear all;              % clears workspace
clc;                    % clears command line

%% Part 1 - Airfoil Interpolation

% Import the data from the 'airfoil_1.txt' document. This data contains 
% the normalized y and z locations of an airfoil surface at 6 locations 
% from trailing to leading edge at the root of the wing. Use the polyfit()
% command to find a 5th order polynomial that interpolates the airfoil.

% Import the data from the 'airfoil_2.txt' document. This data contains 
% the normalized y and z locations of an airfoil surface at 6 locations 
% from trailing to leading edge at the tip of the wing. Use the polyfit() 
% command to find a 5th order polynomial that interpolates the airfoil.

% Using both polynomials from part a and b, use the surf() command to
% plot the wing.

% Airfoil 1
airfoil_1 = load('Lab5_Data\Lab5_Data\airfoil_1.txt');
y_1 = airfoil_1(:,1); % call every row but only first column
z_1 = airfoil_1(:,2); % call every row but only second column

% Airfoil 2
airfoil_2 = load('Lab5_Data\Lab5_Data\airfoil_2.txt');
y_2 = airfoil_2(:,1);
z_2 = airfoil_2(:,2);

% Mesh for Airfoil 1 Y-axis
y_mesh_1 = linspace(y_1(1), y_1(end), 100); % make an array from 1 to whatever max-value of y_1, 100 times in between

% Mesh for Airfoil 2 Y-axis
y_mesh_2 = linspace(y_2(1), y_2(end), 100); % make an array from 1 to whatever max-value of y_2, 100 times in between

% Polyfit & Polyval for Airfoil 1 Z-axis
p_1 = polyfit(y_1,z_1,5); % find coefficients
z_mesh_1 = polyval(p_1,y_mesh_1); % find z height relative to the mesh created

% Polyfit & Polyval for Airfoil 2 Z-axis
p_2 = polyfit(y_2,z_2,5);
z_mesh_2 = polyval(p_2,y_mesh_2);

% Mesh for Airfoil 1 & 2 X-axis
x_mesh = [zeros(100,1), ones(100,1)]; % create an arbitray array to work with when plotting

% Airfoil Plot %

% Set view to near-isometric (solution to figure is slightly off isometric)
v = [2 -1 1]; % coordinates of camera

figure
surf(x_mesh,[y_mesh_1', y_mesh_2'],[z_mesh_1', z_mesh_2']), shading interp % right colors/shade
hold on;
surf(x_mesh,[y_mesh_1', y_mesh_2'],[-z_mesh_1', -z_mesh_2']), shading interp
view(v) % changes view
axis([0 1 0 1 -0.5 0.5]) % set axis boundaries
xlabel('x - Normalized')
ylabel('y - Normalized')
zlabel('z - Normalized')
title('5th Order Polynomial Approximation of Airfoil')

%% Part 2 - Wind Tunnel Data & Regression

% Import the data from the files 'cl_data.txt'  and 'cd_data.txt'. The 
% data is formatted as, first column – angle of attack in degrees, second 
% column – lift/drag coefficient.

% Plot the data from part a and label the figure. Remember the data in the 
% files are data points. It is generally never appropriate to connect
% measured data points with a line.

% Use the polyfit() command to determine the zero angle of attack lift 
% coefficient, 'cl', and the proportionality constant between angle of 
% attack and lift coefficient. In this case, not all the data fits are 
% line. Be sure to limit your regression line to use only points 
% applicable to thin airfoil theory.

% Use the polyfit() command to determine a quadratic fit to the angle of 
% attack versus drag coefficient data.

% Plot both regression lines on the same figure as the data points. Label 
% your figure and be sure to include a legend.

% CL data
cl_data = load('Lab5_Data\Lab5_Data\cl_data.txt');
deg_attack_cl = cl_data(:,1);
rad_attack_cl = deg2rad(deg_attack_cl); % must convert data to radians
lift_cl = cl_data(:,2);

% CD data
cd_data = load('Lab5_Data\Lab5_Data\cd_data.txt');
deg_attack_cd = cd_data(:,1);
rad_attack_cd = deg2rad(deg_attack_cd);
drag_cd = cd_data(:,2);

% Regression Line Code %
% Note that the points used for this line must adhere to the thin airfoil
% theory. This means only columns 1 - 15 as anything above 10 degrees
% violates this theory.

% Regression Line for Lift Data
y_cl_theory = rad_attack_cl(1:15,1);
x_cl_theory = cl_data(1:15,2);

% Regression Line for Lift Data
lift_line = polyfit(y_cl_theory,x_cl_theory,1);

% "lift_line" Equation
x_l = -0.1:0.001:0.35;
y_l = 5.9458*x_l + 0.1533;

% Regression Quadratic Line for Drag Data
drag_line = polyfit(rad_attack_cd,drag_cd,2);

% "drag_line" Equation
x_d = -0.1:0.001:0.35;
y_d = 1.1235*(x_d.^2) - 0.0819*x_l + 0.0151;

% Wind Tunnel Plot w/ Regression Lines %

figure
plot(rad_attack_cl,lift_cl,'o')
hold on
grid on
plot(rad_attack_cd,drag_cd,'o')
plot(x_l,y_l,'Color','b')
plot(x_d,y_d,'Color','r')
xlabel('Angle of Attack (radians)')
ylabel('Coefficient (no units)')
title('Wind Tunnel Data for C-Lift & C-Drag of Airfoil w/ Regression Fits')
legend('C-Lift: Wind Tunnel', 'C-Drag: Wind Tunnel', 'C-Lift: Regression', 'C-Drag: Regression', 'location', 'northwest')
% Davis Sok
% Aero 300 PRELAB5 : 2024.05.01

close all;              % clears all
clear all;              % clears workspace
clc;                    % clears command line

% Define x and y as coordinate points
x = [-5, -2, 4, 5];
y = [4, -1, 2, -5];

% Polyfit command of x and y finding a degree 3 coeffecient
p = polyfit(x,y,3);

% Display results. Each coefficient is assigned to their respective degree
disp('    x^3       x^2       x^1       c')
disp(p)

% Define x and y domain to plot line graph of polyfit
x1 = -5:0.1:5;
y1 = polyval(p,x1);

% Plot figure
figure
plot(x,y,'x','MarkerEdgeColor','r','MarkerSize',12)
hold on
grid on
plot(x1,y1)
xlabel('X Domain')
ylabel('Y Domain')
title('Line of Best Fit via polyfit Command')

% The command polyfit() appears to be most similar to Newton's Divided 
% Difference method in that they are both interpolation methods to find a
% line of best fit w/ a polynomial of a certain degree. As for QR
% Factorization, I unfortunately cannot make any connection. QR
% Factorization seems to relate more to matrices and composing them into a
% lower and upper matrix to solve for A. QR Factorization COULD relate more
% to the Least Squares Method which finds the approximate solution a system
% of equations which would lead to curve fitting.
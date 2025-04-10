%Roshan Jaiswal-Ferri
%Aero 215 Lab 3: 09/28/23

close all;      %Clears all
clear all;      %Clears Workspace
clc;            %Clears Command Window

%variables: 

%Creates an evenly spaces vector starting at 0 going to 2pi,
%and the last number determines how many total values are calculated if
%left blank it will default to 100, replace 2*pi with 360 and sin with sind
%to work with degrees instead.

xValues = linspace(0, 2*pi, 25);

yValues = sin(xValues);

y2Values = cos(xValues);
%Figure creates a popup wondow, which we can put a graph on

figure;
plot(xValues, yValues, '*') %'* showes the plot points with no line. two sashes (--) is a dashed line
xlabel('Angles [rad] from 0 to 2pi');
ylabel('Sin(x)')
title('Graph of Sin(x)')
grid on;
grid minor; %Adds more smaller grid lines

figure;
plot(xValues, y2Values, 'r', 'lineWidth', 3) %the 'r' changes the color to red, 'lineWidth' with the following comma and number changes the width default is one
hold on; %Tells matlab to add another plot on top
plot(xValues, yValues,'g') % Now, both lines will appear
xlabel('Angles [rad] from 0 to 2pi');
ylabel('Cos(x)')
title('Graph of Cos(x)')
grid on;
grid minor; %Adds more smaller grid lines
legend('Sin(x)', 'Cos(x)', 'Location', 'best') %Creates a legend for the graph in order of which plot comes first, the location can also be changed to any cardinal direction as well as 'best' where matlab will autoplace it

figure;
subplot(1, 2, 1); %puts multiple plots on one figure/window with the numbers changing the position
plot(xValues, yValues, '*') %'* showes the plot points with no line. two sashes (--) is a dashed line
xlabel('Angles [rad] from 0 to 2pi');
ylabel('Sin(x)')
title('Graph of Sin(x)')
grid on;
grid minor; %Adds more smaller grid lines
subplot(1, 2, 2)
plot(xValues, y2Values, 'r', 'lineWidth', 3) %the 'r' changes the color to red, 'lineWidth' with the following comma and number changes the width default is one
hold on; %Tells matlab to add another plot on top
plot(xValues, yValues,'g') % Now, both lines will appear
xlabel('Angles [rad] from 0 to 2pi');
ylabel('Cos(x)')
title('Graph of Cos(x)')
grid on;
grid minor; %Adds more smaller grid lines











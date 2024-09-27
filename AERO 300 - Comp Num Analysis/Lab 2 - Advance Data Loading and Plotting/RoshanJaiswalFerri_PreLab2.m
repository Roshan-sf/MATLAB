%Roshan Jaiswal-Ferri
%Section - 03
%Aero 300 Pre Lab 2 - Advance Data Loading and Plotting: 4/9/24

%%

close all;      %Clears all
clear all;      %Clears Workspace
clc;            %Clears Command Window

%% PART 1: plot()

o = -2*pi; %Setting bounds to variables
p = 2*pi;

theta = linspace(o,p,130); %Creating vector with bounds and stepping

y = pi*sin(theta/2); %Example Function

g = theta/2; %Other Example function

figure; %Creating a figure with overlayed functions using plot
plot(theta, y)
hold on
plot(theta, g) %using the plot command to create two overlaying lines on single figure
grid on;
title('plot(): Graph of πSin(θ/2) & θ/2')

%% PART 2: contour()

figure;
x = linspace(-2*pi,2*pi);
y = linspace(0,4*pi);
[X,Y] = meshgrid(x,y);
Z = sin(Y)+cos(X);
contour(X,Y,Z) %the contour function creates a topographic map based on x y coordinates with height z
title('contour()')

%% PART 3: surf()

figure;
surf(X,Y,Z) %using the same variables used for the contour command, surf creates an actual 3d graph of x y z data
title('surf()')

%% Part 4: streamline()

load wind %Using built in wind vector data for this example
[startX,startY,startZ] = meshgrid(80,20:10:50,0:5:15); %Setting start points
verts = stream3(x,y,z,u,v,w,startX,startY,startZ);
figure
lineobj = streamline(verts);
view(3)
title('streamline()')

%figure
%streamline(X,Y,Z,U,V,W,startX,startY,startZ) %Creates lines following 2D or 3D vector data

%% Part 5: quiver()

load('wind','x','y','u','v')
X = x(11:22,11:22,1);
Y = y(11:22,11:22,1);
U = u(11:22,11:22,1);
V = v(11:22,11:22,1);

figure
quiver(X,Y,U,V) %Creates a 2D vector plot of given data (Same wind data as earlier)
title('quiver()')



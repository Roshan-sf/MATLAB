clc
clear
close all




psi = 0;
theta = 0;
phi = 0;


xrot = [1,0,0;0,cos(psi),sin(psi);0,-sin(psi),cos(psi)];
yrot = [cos(theta) 0 -sin(theta);0 1 0;sin(theta) 0 cos(theta)];
zrot = [cos(phi) sin(phi) 0;-sin(phi) cos(phi) 0;0 0 1];

%% Part 1

psi = pi;
Cenuac0 = [1,0,0;0,cos(psi),sin(psi);0,-sin(psi),cos(psi)];
disp(Cenuac0)

%% Part 2.1 
psi = (20*pi)/180;
theta = (-5*pi)/180;
phi = (10*pi)/180;




Cac0acf = [1,0,0;0,cos(psi),sin(psi);0,-sin(psi),cos(psi)]*[cos(theta) 0 -sin(theta);0 1 0;sin(theta) 0 cos(theta)]*[cos(phi) sin(phi) 0;-sin(phi) cos(phi) 0;0 0 1];
disp(Cac0acf)

%% Part 2.2
psi = (20*pi)/180;
theta = (-5*pi)/180;
phi = (10*pi)/180;
xrot = [1,0,0;0,cos(psi),sin(psi);0,-sin(psi),cos(psi)];
yrot = [cos(theta) 0 -sin(theta);0 1 0;sin(theta) 0 cos(theta)];
zrot = [cos(phi) sin(phi) 0;-sin(phi) cos(phi) 0;0 0 1];

Flipped = zrot*yrot*xrot;
Cacfenu = Cenuac0*Flipped;

%2.3

Northenu = [0; 1 ;0];
Northofac = Cacfenu*Northenu;







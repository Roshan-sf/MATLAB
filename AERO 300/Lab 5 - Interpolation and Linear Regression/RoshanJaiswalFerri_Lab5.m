%Roshan Jaiswal-Ferri
%Section - 03
%Aero 300 Lab 5 - Interpolation and Linear Regression: 5/2/24

%% Workspace Prep

format long     %Allows for more accurate decimals
close all;      %Clears all
clear all;      %Clears Workspace
clc;            %Clears Command Window

%% PART 1: Plot the Wing — Interpolation

%Airfoil 1
af1 = load("Lab5_Data\Lab5_Data\airfoil_1.txt");
y1 = af1(:,1); %every row in column 1
z1 = af1(:,2); %every row in column 2

%Airfoil 2
af2 = load("Lab5_Data\Lab5_Data\airfoil_2.txt");
y2 = af2(:,1);
z2 = af2(:,2);

%Spacing Vector for Y-axis
yVec1 = linspace(y1(1), y1(end), 100); 
yVec2 = linspace(y2(1), y2(end), 100); 

%Creating polynomials for z:
p1 = polyfit(y1,z1,5); %coefficients of degree 5 poly af1
p2 = polyfit(y2,z2,5); 

zVal1 = polyval(p1,yVec1); %z height
zVal2 = polyval(p2,yVec2);

x = [zeros(100,1), ones(100,1)]; %Vector to plot on

%Plotting

figure('Name','Wing Plot')
surf(x,[yVec1', yVec2'],[zVal1', zVal2']), shading interp %proper shading
hold on;
surf(x,[yVec1', yVec2'],[-zVal1', -zVal2']), shading interp
xlabel('X')
ylabel('Y')
zlabel('Z')
title('Approx. of Airfoil with 5th Degree Polynomial')
linestyle = 'none';
axis([0 1 0 1 -0.5 0.5]) % set axis boundaries
view([2 -1 1]) %correct view

%% PART 2: Wind Tunnel Data — Regression 

%CL data
cl = load('Lab5_Data\Lab5_Data\cl_data.txt');
dAcL = cl(:,1); %degree angle of attack CL 
rAcL = deg2rad(dAcL); %radian angle of attack CL
cll = cl(:,2); %CL lift

%CD data
cd = load('Lab5_Data\Lab5_Data\cd_data.txt');
dAcD = cd(:,1); %degree angle of attack CD
rAcD = deg2rad(dAcD); %radian angle of attack CD
cdd = cd(:,2); %CD drag

%%%%Code used for finding regression for lift and drag%%%%

%Lift regression
ycL = rAcL(1:15,1);
xcL = cl(1:15,2);
lift = polyfit(ycL,xcL,1); %Regression line

%Drag regression
drag = polyfit(rAcD,cdd,2);

%%%%

%lift eq
x_l = -0.1:0.001:0.35;
y_l = 5.9458*x_l + 0.1533;

%drag eq
x_d = -0.1:0.001:0.35;
y_d = 1.1235*(x_d.^2) - 0.0819*x_l + 0.0151;

%Plotting
figure('Name','Wind Tunnel Regression for Coeff. of Drag and Lift')
title('Wind Tunnel Regression for Coeff. of Drag and Lift')
plot(rAcD,cdd,'o','color','r')
hold on
grid on
plot(rAcL,cll,'o','color','g')
plot(x_l,y_l,'color','g')
plot(x_d,y_d,'color','r')
xlabel('Angle of Attack (rad)')
ylabel('Coeff. (unitless)')
legend('Coeff. Drag', 'Coeff. Lift', 'Coeff. Lift Regression', 'Coeff. Drag Regression', 'location', 'best')










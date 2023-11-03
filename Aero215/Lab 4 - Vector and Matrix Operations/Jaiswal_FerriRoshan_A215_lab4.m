%Roshan Jaiswal-Ferri
%Aero 215 Lab 4: 10/31/23

%%

close all;      %Clears all
clear all;      %Clears Workspace
clc;            %Clears Command Window

%% Part 1 Defining Variables

H = [3,2,1];
N = [9,8,7];

M = [H;N]; %creates a matrix ; means go down a row , is next column

V = [1,2,3,4,5];
t = [0.2,0.4,0.6,0.9,1.1];

%% Part 2: Cross Product

R = [(M(1,2)*M(2,3))-(M(1,3)*M(2,2)), -((M(1,1)*M(2,3))-(M(1,3)*M(2,1))), (M(1,1)*M(2,2))-(M(1,2)*M(2,1))]; %Do a cross product

disp(num2str(R)) %Display result
disp([num2str(cross(H,N))]) %Check if correct (it is)

%% Part 3: Vector Graphs

V2 = V.^2;
V3 = V.^3;

figure;
plot(t,V, 'r');
hold on; 
plot(t, V2,'g')
plot(t, V3,'y')
xlabel('Time (s)');
ylabel('Velocity')
title('Velocity vs Time')
grid on;
legend('V', 'V^2', 'V^3', 'Location', 'best')



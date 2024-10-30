%Roshan Jaiswal-Ferri
%Section - 01
%Aero 320 - 10/29/24

%% Workspace Prep

format long     %Allows for more accurate decimals
close all;      %Clears all
clear all;      %Clears Workspace
clc;            %Clears Command Window

%% PART 1:

%J = A really big matrix

initialAngVec = [0; 0; 0.005];

initialEA = [0; 0; -deg2rad(45)];

initialTorque = [0; 0; 490.50];

J 

%%

initialStates = [initialAngVec; initialEA];
tspan = [0 45];

[t, states] = ode45(@(t,states) evolution(t, initialStates, J, initialTorque), tspan, initialStates, []);

%%


function dy = evolution(t,y,J,T)

    %Initialize column - ang vel, eul ang
    dy(1:6,1) = zeros(6,1);

    %Extract states
    w = y(1:3,1);
    phi = y(4,1);
    theta = y(5,1);
    psi = y(6,1);
    
    %propogate ang vel
    %Euler's equation of motion (torques)
    dx(1:3,1) = inv(J)*(T-vect2cross(w)*J*w);

    %Euler Rates
    mtx = [cos(theta), sin(phi)*sin(theta), cos(phi)*sin(theta);...
        0, cos(phi)*cos(theta), -sin(phi)*cos(theta);...
        0, sin(phi), cos(phi)];
    dy(4:6,1) = 1/cos(theta)*mtx*w;

end

function [out] = vect2cross(v)
    out = [0, -v(3), v(2); v(3), 0, -v(1); -v(2), v(1), 0];
end


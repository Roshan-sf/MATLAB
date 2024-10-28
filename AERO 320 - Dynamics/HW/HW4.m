%Roshan Jaiswal-Ferri
%Section - 01
%Aero 320 HW 4 - Problem 2: 10/21/24

%% Workspace Prep

format long     %Allows for more accurate decimals
close all;      %Clears all
clear all;      %Clears Workspace
clc;            %Clears Command Window

%% PART 1:

w21 = [1; -2; .5]; %angular vel
q0 = [0; 0; 0; 1]; %initial quaternion
eulA = [0; 0; 0]; %initial euler angles of 0

time = [0, 10]; %sec

[t1, q] = ode45(@quaternion, time, q0, [], w21); %
[t2, e] = ode45(@eula, time, eulA, [], w21);

figure('Name','Evolution of Euler Angles')
plot(t2, e(:,1), t2, e(:,2), t2, e(:,3));
grid on
xlabel('Time (s)');
ylabel('Euler Angles (deg)');
legend('phi', 'theta', 'psi');
title('Evolution of Euler Angles');

figure('Name','Evolution of Quaternions')
plot(t1, q(:,1), t1, q(:,2), t1, q(:,3), t1, q(:,4));
grid on
xlabel('Time (s)');
ylabel('Quaternion Values');
legend('E1', 'E2', 'E3', 'eta');
title('Evolution of Quaternions');


%% Functions:

function [dq] = quaternion(t,x,o)
    ep = x(1:3,1);
    n = x(4);

    dq(1:3,1) = .5*(n*eye(3)+vcross(ep))*o;
    dq(4,1) = -.5*ep'*o;
end

function [de] = eula(t,x,o)
    ph = x(1);
    th = x(2);

    de(1:3,1) = 1/cos(th)*[cos(th), sin(ph)*sin(th), cos(ph)*sin(th);...
        0, cos(ph)*cos(th), -sin(ph)*cos(th);...
        0, sin(ph), cos(ph)]*o;
end

function [out] = vcross(v)
    out = [0, -v(3), v(2); v(3), 0, -v(1); -v(2), v(1), 0];
end
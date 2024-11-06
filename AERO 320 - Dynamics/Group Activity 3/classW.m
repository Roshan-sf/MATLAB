%Roshan Jaiswal-Ferri
%Section - 01
%Aero 320 HW 5 - 10/29/24

%% Workspace Prep

format long;  % Allows for more accurate decimals
close all;    % Clears all
clear all;    % Clears Workspace
clc;          % Clears Command Window

%% PART 1: Problem 1 Part A

wb0 = [0.05; 0.01; 3.5]; %angular vel rad/s
q0 = [0; 0; 0; 1]; %initial quaternion
eulA = [0; 0; 0]; %initial euler angles of 0
Torque = [0; 0; 0]; %initial torque in Nm
time = [0, 20]; %sec

J = [45, 0, 0; 0, 45, 0; 0, 0, 70];

%-10*sin(deg2rad(.5))
Torque = [0; 0; 0];  % Applied torque in Nm

[t3, w] = ode45(@(t, states) angularV(t, states, J, Torque), time, wb0);
w_interp = @(t) interp1(t3, w, t);

[t1, q] = ode45(@(t, y) quaternion(t, y, w_interp(t).'), time, q0);
[t2, e] = ode45(@(t, y) eula(t, y, w_interp(t).'), time, eulA);

figure('Name', 'Roll Pitch & Yaw');
% subplot(3, 1, 1);  % 3 rows, 1 column, 1st subplot
plot(t2, e(:,1), t2, e(:,2), t2, e(:,3));
grid on;
xlabel('Time (s)');
ylabel('Euler Angles (deg)');
legend('\phi', '\theta', '\psi');
title('Evolution of Euler Angles');

% subplot(3, 1, 2);  % 3 rows, 1 column, 2nd subplot
figure('Name', 'Quaternion Propogation');
plot(t1, q(:,1), t1, q(:,2), t1, q(:,3), t1, q(:,4));
grid on;
xlabel('Time (s)');
ylabel('Quaternion Values');
legend('E1', 'E2', 'E3', '\eta');
title('Evolution of Quaternions');

figure('Name', 'Angular Velocity Over Time');
%subplot(3, 1, 3);  % 3 rows, 1 column, 3rd subplot
plot(t3, w(:,1), t3, w(:,2), t3, w(:,3));
grid on;
xlabel('Time (s)');
ylabel('Angular Velocity (rad/s)');
legend('\omega_x', '\omega_y', '\omega_z');
title('Angular Velocity Over Time');

%% PART 2: Problem 1 part b (with torque)








%% Functions:

function dw = angularV(t,w,J,T)

    %Euler's equation of motion (torques)
    dw(1:3,1) = inv(J)*(T-vcross(w)*J*w);

end

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

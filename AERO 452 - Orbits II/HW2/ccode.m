%% Problem 2
clc
clear
close all
r0 = [0, 0, 0]';
T = 24*60*60;
n = 2*pi/T; %rad/s
t = 2 * 60 * 60;
phi_rr = [4-3*cos(n*t), 0, 0;
    6*(sin(n*t)-n*t), 1, 0;
    0, 0, cos(n*t)];
phi_rv = [1/n*sin(n*t), 2/n*(1-cos(n*t)), 0; ...
    2/n*(cos(n*t)-1), 1/n*(4*sin(n*t)-3*n*t), 0;
    0, 0, 1/n *sin(n*t)];
phi_vr = [3*n*sin(n*t), 0, 0;
    6*n*(cos(n*t)-1), 0, 0;
    0, 0, -n*sin(n*t)];
phi_vv = [cos(n*t), 2*sin(n*t), 0;
    -2*sin(n*t), 4*cos(n*t)-3, 0;
    0, 0, cos(n*t)];
r2 = [-10, 10, 0]'.*1000;
v0 = inv(phi_rv)*(r2-phi_rr*r0);
v2 = phi_vr*r0 + phi_vv *v0;


t6 = 6*60*60;
[deltav0_plus, deltav_final_minus] = CW_matrix(n, t6, r2);
deltav_total = norm(deltav0_plus - v2) + norm(deltav_final_minus);
disp(['Problem 2']);
disp(['Delta-v: ', num2str(deltav_total), ' m/s']);

function [deltav0_plus, deltav_final_minus] = CW_matrix(n, t, r0)
phi_rr = [4-3*cos(n*t), 0, 0;
    6*(sin(n*t)-n*t), 1, 0;
    0, 0, cos(n*t)];
phi_rv = [1/n*sin(n*t), 2/n*(1-cos(n*t)), 0; ...
    2/n*(cos(n*t)-1), 1/n*(4*sin(n*t)-3*n*t), 0;
    0, 0, 1/n *sin(n*t)];
phi_vr = [3*n*sin(n*t), 0, 0;
    6*n*(cos(n*t)-1), 0, 0;
    0, 0, -n*sin(n*t)];
phi_vv = [cos(n*t), 2*sin(n*t), 0;
    -2*sin(n*t), 4*cos(n*t)-3, 0;
    0, 0, cos(n*t)];
deltav0_plus = inv(phi_rv)*-phi_rr*r0;
deltav_final_minus = phi_vr * r0 + phi_vv * deltav0_plus;
end
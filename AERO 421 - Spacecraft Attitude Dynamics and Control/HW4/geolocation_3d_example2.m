clear
close all
clc

c = [0;0;0];
Re = 6378;

%% Simple Example - Looking back Nadir
% Define nominal pointing direction.  Looking directly back nadir
p_dir = [-1; 0; 0];
p_hat = p_dir/norm(p_dir);

% Location of spacecraft
r_0 = [7378;0;0];

% Solve for nominal target location
[T_0, d] = sphere_line_intersect(c, Re, p_hat, r_0)

%% More complex example - looking back at some angle

% Define nominal pointing direction.  Mostly directly back nadir
p_dir_e = [-1; .01; -.02];
p_hat_e = p_dir_e/norm(p_dir_e);

% Location of spacecraft
r_0 = [7378;0;0];

% Solve for nominal target location
[T_e, d] = sphere_line_intersect(c, Re, p_hat_e, r_0)


ground_err = norm(T_0 - T_e)
theta_err = acos(p_hat'*p_hat_e)*180/pi
% figure
% [ex, ey, ez] = ellipsoid(0, 0, 0, Re, Re, Re, 30);
% sEarth = surf(ex, ey, flip(ez), 'EdgeColor','none', 'FaceAlpha', .7);
% r = arrow([0 r_0(1)], [0 r_0(2)], [0 r_0(3)], 'facecolor', 'red');
% axis vis3d
% hold on;


% Number of Monte Carlo iterations
n = 1000;

% Preallocate error and target matricies
e = zeros(3,n);
e_norm = zeros(1,n);
T = zeros(3,n);

% Set the random number generator
rng('default')
s = rng;

% Begin Monte Carlo Loop
for i = 1:n
    % calc all sources of error for r, d_r_pos, d_r_sen, d_r_vel
    d_r_pos = 10*randn([3, 1]);
    r = r_0 + d_r_pos;
    
    % calc the errors for target location (d_lat, d_long, d_alt, d_JD)
    %d_T(:,i) = lla2eci(lat_0 + d_lat, long_0 + d_long, alt_0 + d_alt, d_JD);
    T(:,i) = T_0;% + d_T;
    
    % calc all erros for the pointing vector - sensor mounting angle
    % d_p = ...
    % p_hat = p_hat_0 + d_p
    
    % Solve for the location where p intersects Earth
    T_e = sphere_line_intersect(c, Re, p_hat, r);

    % Find Ground Error
    e(:,i) = T(:,i) - T_e;
    e_norm(i) = norm(T(:,i) - T_e);
end

% Gather Point Error Statistics
std(e_norm);
mean(e_norm); %km

e_s = sort(abs(e_norm));
e90 = e_s(floor(n*0.9));

%% Make Plots
figure
plot(e_norm);
xlabel('trial')
ylabel('distance error (km)')
title(horzcat('random error of ', num2str(n), ' trials'))

figure
plot(sort(abs(e_norm)))
grid on
xlabel('trial')
ylabel('distance error (km)')
title(horzcat('sorted absolute random error of ',num2str(n),' trials'))

figure
plot3(e(1,:)+T(1,1), e(2,:), e(3,:),'.')
grid on
title(horzcat('Scatter plot of ground error for ', num2str(n), ' iterations.'))
zlabel('error in z (km)')
ylabel('error in y (km)')
view([-90,0])


function [T_0, d] = sphere_line_intersect(c, Re, p_hat, r_0)
    r = r_0 - c;
    
    % Compute quadratic terms
    
    discriminant = dot(r, p_hat)^2 - (dot(r,r) - Re^2);
    
    if discriminant < 0
        warning('Discriminant less than zero, result will be imag')
        d = 0;
        T_0 = [0,0,0];
    else
        d1(1) = dot(-r,p_hat)+sqrt(discriminant);
        d1(2) = dot(-r,p_hat)-sqrt(discriminant);
        d = min(d1);

        T_0 = r + d*p_hat;
    end
end

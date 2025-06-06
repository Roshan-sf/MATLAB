%% Roshan Jaiswal-Ferri
%Section - 01
%Aero 421 HW4: 5/9/25

%% Workspace Prep

%warning off
format long     %Allows for more accurate decimals
close all;      %Clears all
clear all;      %Clears Workspace
clc;            %Clears Command Window

%% Variables

Re = 6378; %km
mu = 398600;
JD = 2458981.1666667;
c = [0; 0; 0];
n = 1000;

r0 = [5980.8297; -1282.3184; 4125.8019];  % km in ECI
v0 = [1.540887; 7.186813; 0];  % km/s

target.lat = 35.3;       % degrees
target.lon = -120.8;     % degrees
target.alt = 0.2;        % km

target_eci = lla2eci(target.lat, target.lon, target.alt, JD); %T_0
T_0 = target_eci;

p_0 = T_0 - r0;
phat_0 = p_0/norm(p_0);


%-------------------------------ERRORS-------------------------------------
% Local Sidereal Time
error.lsdErr = 0.01; % seconds
error.JDError = 0.01/86400; % now in days which is jd or somethin idk

% Target Location
error.LatErr = deg2rad(1e-4); % degrees to read
error.LonErr = deg2rad(1e-4); % degrees to rad
error.altErr = 10/1000; % kilometers

% Spacecraft Knowledge
error.rErr = 3/1000; % meters (in-track, cross-track, radius)
error.vErr = 0.002/1000; % m/s (in-track, cross-track)
error.vErr_r = 0.007/1000; % m/s (radial)

% Sensor Mounting Errors
error.smErr = 0.01/1000; % meters
error.smaErr = deg2rad(1e-4); % degrees to rad

%% Problem 1

[g, g_norm, T_g,e] = geo_monte_sim(r0, v0, phat_0, T_0, JD, error, target, Re, n, c);

disp('Part 1:')
disp(['g vector in meters: ', num2str(g')])
disp(['gnorm: ', num2str(g_norm')])
disp(['Tg: ', num2str(T_g')])
disp(' ')

%% Problem 2

altVec = linspace(1000,9000,20);
g_errors = zeros(length(altVec), 1);
slant_ranges = zeros(length(altVec), 1);

rhat = r0/norm(r0);
vhat = v0/norm(v0);

for i = 1:length(altVec)
    rmag = Re + altVec(i);
    vmag = sqrt(mu/rmag);
    
    r02 = rmag*rhat;
    v02 = vmag*vhat;

    p_02 = T_0 - r02; %update pointing vector
    phat_02 = p_02/norm(p_02);

    [~, loc_error, ~, mcData] = geo_monte_sim(r02, v02, phat_02, T_0, JD, error, target, Re, n, c);
   
    g_errors(i) = loc_error;
    slant_ranges(i) = norm(T_0 - r02);
    mc_results{i} = mcData;
end

figure;
plot(altVec, g_errors, '-o');
xlabel('Altitude (km)');
ylabel('90% Geolocation Error (m)');
title('Geolocation Error vs. Altitude');

figure;
plot(slant_ranges, g_errors, '-o');
xlabel('Slant Range (km)');
ylabel('90% Geolocation Error (m)');
title('Geolocation Error vs. Slant Range');

% Save data for external analysis
save('hw4_q2_results.mat', ...
     'altVec', ...
     'g_errors', ...
     'slant_ranges', ...
     'rhat', ...
     'vhat', ...
     'Re', ...
     'mu', ...
     'JD', ...
     'target', ...
     'T_0', ...
     'error', ...
     'n', ...
     'c', ...
     'mc_results');

%% Functions


function [g, g_norm, T_g, monteCarloData] = geo_monte_sim(r_0, v_0, phat_0, T_0, JD, error, target, Re, n, c)
    % Target location
    lat_0 = target.lat;
    long_0 = target.lon;
    alt_0 = target.alt;

    % Local Sidereal Time
    lsdErr = error.lsdErr;
    JDErr = error.JDError;

    % Target Location error
    LatErr = error.LatErr;
    LonErr = error.LonErr;
    altErr = error.altErr;

    % Spacecraft Knowledge
    rErr = error.rErr;
    vErr = error.vErr;
    vErr_r = error.vErr_r;

    % Sensor Mounting Errors
    smErr = error.smErr;
    smaErr = error.smaErr;

    rng('default')

    for i = 1:n
        % Apply Gaussian noise to position and velocity
        d_r_pos = error.rErr * 1000 * randn(3,1);  % Apply in meters, convert to km
        d_r_vel = error.vErr * 1000 * lsdErr * randn(3,1);
        d_r_radial = (error.vErr_r * 1000 * lsdErr * randn()) * (r_0 / norm(r_0));

        r = r_0 + d_r_pos + d_r_vel + d_r_radial;

        % Perturbed target location
        tlat = lat_0 + LatErr * randn();
        tlon = long_0 + LonErr * randn();
        talt = alt_0 + altErr * randn();
        tJD = JD + JDErr * randn();

        T(:,i) = lla2eci(tlat, tlon, talt, tJD);

        % Local sensor frame based on current r
        phatz = phat_0;
        r_unit = r / norm(r);
        phaty = cross(r_unit, phatz); phaty = phaty / norm(phaty);
        phatx = cross(phaty, phatz); phatx = phatx / norm(phatx);

        % Apply small-angle rotation (sensor mounting error)
        phi = smaErr * randn();  % angular deviation
        theta = 2 * pi * rand(); % random direction

        d_p = phi * (cos(theta) * phatx + sin(theta) * phaty);
        d_p = d_p + error.smErr * 1000 * randn(3,1);

        p_hat = (phat_0 + d_p) / norm(phat_0 + d_p);

        % Line-of-sight intersection
        T_e(:,i) = sphere_line_intersect(c, Re, p_hat, r);

        % Error vector
        e(:,i) = T(:,i) - T_e(:,i);
        e_norm(i) = norm(e(:,i));
    end

    z90 = 1.645;
    g = mean(e, 2) + z90 * std(e, 0, 2);
    g_norm = mean(e_norm) + z90 * std(e_norm);
    T_g = mean(T_e, 2) + z90 * std(T_e, 0, 2);

    % Save internal data
    monteCarloData.e = e;
    monteCarloData.e_norm = e_norm;
    monteCarloData.T_e = T_e;
    monteCarloData.T = T;
end

function [T_0, d] = sphere_line_intersect(c, Re, p_hat, r_0)
    r = r_0 - c;
       
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


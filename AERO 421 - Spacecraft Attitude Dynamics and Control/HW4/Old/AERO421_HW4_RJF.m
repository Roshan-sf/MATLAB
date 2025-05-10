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

    [~, loc_error] = geo_monte_sim(r02, v02, phat_02, T_0, JD, error, target, Re, n, c);
   
    g_errors(i) = loc_error;
    slant_ranges(i) = norm(T_0 - r02);
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

%% Functions

function [g, g_norm, T_g,e] = geo_monte_sim(r_0, v_0, phat_0, T_0, JD, error, target, Re, n, c)
    % Target location
    lat_0 = target.lat;
    long_0 = target.lon;
    alt_0 = target.alt;

    % Local Sidereal Time
    lsdErr = error.lsdErr; % seconds
    JDErr = error.JDError; %now in days which is jd or somethin idk
    
    % Target Location error
    LatErr = error.LatErr; % degrees to read
    LonErr = error.LonErr; % degrees to rad
    altErr = error.altErr; % meters
    
    % Spacecraft Knowledge
    rErr = error.rErr; % meters (in-track, cross-track, radius)
    vErr = error.vErr; % m/s (in-track, cross-track)
    vErr_r = error.vErr_r; % m/s (radial)
    
    % Sensor Mounting Errors
    smErr = error.smErr; % meters
    smaErr = error.smaErr; % degrees to rad

    % Set the random number generator
    rng('default')
    s = rng;
    
    % Begin Monte Carlo Loop
    for i = 1:n
        % calc all sources of error for r, d_r_pos, d_r_sen, d_r_vel
        d_r_pos = rErr*randn([3, 1]);
        d_lat = LatErr*rand([1, 1]);
        d_long = LonErr*rand([1, 1]);
        d_alt = (altErr)*rand([1, 1]);
        d_JD = JDErr*rand([1, 1]);
    
        d_r_pos = d_r_pos +  + vErr * lsdErr * rand(3,1);
        d_r_pos = d_r_pos + (vErr_r * lsdErr * rand()) * (r_0 / norm(r_0));

        r = r_0 + d_r_pos;
        
        % calc the errors for target location (d_lat, d_long, d_alt, d_JD)
        tlat = lat_0 + d_lat;
        tlon = long_0 + d_long;
        talt = alt_0 + d_alt;
        tJD = JD + d_JD;

        d_T(:,i) = lla2eci(tlat, tlon, talt, tJD);
        T(:,i) = d_T(:,i); %no need to add T_0 b/c it is already considered with inout args
        
        % calc all erros for the pointing vector - sensor mounting angle

        phatz = phat_0;
        phaty = cross((r_0/norm(r_0)),phat_0);
        phatx = cross(phaty,phatz);

        phi = randn([1, 1]);
        theta = smaErr*rand([1, 1]);
        d_p = (phi*cos(theta)*phatx) + (phi*sin(theta)*phaty);
        d_p = d_p + smErr * rand(3,1);
    
        p_hat = (phat_0 + d_p)/norm((phat_0 + d_p));
        
        % Solve for the location where p intersects Earth
        T_e(:,i) = sphere_line_intersect(c, Re, p_hat, r);
    
        % Find Ground Error
        e(:,i) = T(:,i) - T_e(:,i);
        e_norm(i) = norm(T(:,i) - T_e(:,i));
    end

    z90 = 1.645;

    xbarx = mean(e(1,:));
    xbary = mean(e(2,:));
    xbarz = mean(e(3,:));
    xbar = [xbarx; xbary; xbarz];

    sigmax = std(e(1,:));
    sigmay = std(e(2,:));
    sigmaz = std(e(3,:));
    sigma = [sigmax; sigmay; sigmaz];

    x90 = xbar + z90*sigma;
    g = x90;
    
    xbar2 = mean(e_norm(1,:));

    sigma2 = std(e_norm(1,:));

    x290 = xbar2 + z90*sigma2;
    g_norm = x290;

    xbarx = mean(T_e(1,:));
    xbary = mean(T_e(2,:));
    xbarz = mean(T_e(3,:));
    xbar3 = [xbarx; xbary; xbarz];

    sigmax = std(T_e(1,:));
    sigmay = std(T_e(2,:));
    sigmaz = std(T_e(3,:));
    sigma3 = [sigmax; sigmay; sigmaz];

    x390 = xbar3 + z90*sigma3;
    T_g = x390;

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


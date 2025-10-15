%% Roshan Jaiswal-Ferri
%Aero 452 Homework 2: 10/8/25

%% Workspace Prep

format long     %Allows for more accurate decimals
close all;      %Clears all
clear all;      %Clears Workspace
clc;            %Clears Command Window

%% PART 1: CW Motion

dr0 = [0; 0; 0]; % m/s
dv0 = [1; -1; 1]; % m/s

% n = (2*pi)/T
% t = (1/4)*T
% nt = ( (2*pi)/T ) * ( (1/4)*T ) = pi/2
nt = pi/2;

[~, ~, phi_vr, phi_vv] = SolveCW(nt);

dv = (phi_vr * dr0) + (phi_vv * dv0);
disp(num2str(dv'))
disp(['Q1: ', num2str(norm(dv)), ' m/s'])
disp(' ')

disp('Discussion for Q1:')
fprintf(['This answer makes sense because the motion is periodic, at 1/4 of \n' ...
    ' the period the x and y components combined and the z velocity \n' ...
    '(out of plane) reaches zero as it oscilates and changes direction '])
disp(' ')
disp(' ')

%% Curtis 7.15

r0 = [0; 0; 0]; %origin
r2 = [-10; 10; 0]*1000; %current pos km
tf2 = 2*3600; %time of flight after being hit
tf6 = 6*3600; %flight time in seconds

T = 24*3600; %assuming 24hr period b/c geo orbit
n = (2*pi)/T;

[phirr, phirv, phivr, phivv] = SolveCW(n,tf2);

Vinitial = (phirv) \ (r2-phirr*r0); %right after being hit v0 plus
Vdrift = (phivr * r0) + (phivv * Vinitial); %velocity 2hrs after

[~, ~, ~, ~, dvplus, dvminus] = SolveCW(n,tf6,r2);

dv_total = norm(dvplus - Vdrift) + norm(dvminus);  % m/s
disp(['Q2: ', num2str(dv_total), ' m/s'])
disp(' ')

disp('Discussion for Q2:')
fprintf(['The deltav is small because GEO’s mean motion is low and the \n' ...
    'repositioning occurs over a long period of time.\n\n'])

%% Two-Impulse Target/Chaser

Re = 6378;
r = 300 + Re;
mu = 398600;

r0c = [-1; 0; 0];
tf = 30*60;
P = 2*pi*sqrt((r^3)/mu);
n = (2*pi)/P;

[phi_rr, phi_rv, phi_vr, phi_vv, dvminus, dvplus] = SolveCW(n,tf,r0c);

dv_total = norm(dvplus) + norm(dvminus);  % m/s
disp(['Q3: ', num2str(dv_total*1000), ' m/s'])
disp(' ')

disp('Discussion for Q3:')
fprintf(['Both the dvplus and dvminus vectors have a zero out of plane \n' ...
    'component, this allows the overall total deltav to remain small \n' ...
    'for this rendeszvous.\n\n'])

%% Coplanar Relative Motion

ra = [8000; 0; 0]; %km
rb = [7000; 0; 0];

va = [0; sqrt(mu/norm(ra)); 0];
vb = [0; sqrt(mu/norm(rb)); 0];

rho = rb-ra;
dx = -1000;

h = cross(ra,va);
omega = h/norm(ra)^2;

rhodot = vb - va - cross(omega, rho);
Q = ECI2LVLH(ra,va);

vLVLH = Q*rhodot;

P = 2*pi*sqrt((norm(ra)^3)/mu);
n = (2*pi)/P;

rv = [0, norm(-3/2*n*dx), 0];

disp('Q4 Describing Relative Motion:')
fprintf(['Because B is lower than A it has a higher angular rate, at perigee \n' ...
    'b is directly below a, there is only a relative distance in tangential \n' ...
    'velocity, since in this frame there is no out of plane or radial change \n\n'])

disp('Q4 Relateive Velocity:')
disp(['LVLH: ', num2str(vLVLH')])
disp(['Relative Motion: ', num2str(rv)])

%% Functions:

function [Phi_rr, Phi_rv, Phi_vr, Phi_vv, dv0plus, dv0minus] = SolveCW(n, t, r0)
%Solves the CW equations in Matrix form
%INPUTS: SolveCW(n,t,r), or (n,t), or (n*t)
%INPUTS: SolveCW(mean motion (rad/s), elapsed time (scaler number), rvec)
%INPUTS: SolveCW(n*t, r0)
%OUTPUTS: [Phi_rr, Phi_rv, Phi_vr, Phi_vv, dv0plus, dv0minus]

    if nargin == 1
        nt = n;
    elseif nargin == 2 || nargin == 3
        nt = n*t;
    end

    s = sin(nt);
    c = cos(nt);

    Phi_rr = [ 4-3*c,           0,   0;
               6*(s-nt),        1,   0;
               0,               0,   c ];

    Phi_rv = [ (1/n)*s,         (2/n)*(1-c), 0;
               (2/n)*(c-1),     (1/n)*(4*s-3*nt), 0;
               0,               0,       (1/n)*s ];

    Phi_vr = [ 3*n*s,           0,   0;
               6*n*(c-1),       0,   0;
               0,               0,  -n*s ];

    Phi_vv = [ c,               2*s, 0;
              -2*s,        4*c-3,   0;
               0,               0,   c ];
    
    if nargin == 3
        dv0plus = (Phi_rv) \ -Phi_rr*r0;
        dv0minus = Phi_vr * r0 + Phi_vv * dv0plus;
    else
        dv0plus = NaN;
        dv0minus = NaN; 
    end

end

